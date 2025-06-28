#! /bin/bash
set -eu -o pipefail

DIR=""
RERUN=false
BINNING_BED=""
BED_GREP_PATTERN="^hor_"
REF=""
READS=""
PAF=""
JOIN=5000000
MARGIN=500000
MIN_READ_LEN=50000
MIN_AVG_IDENT=0.9
MIN_SPAN=10000
MODE="minimap2"
MINIMAP2="mm2-ivh"
MINIASM="miniasm"
HIFIASM="hifiasm"
AVA_OPTS="-xivh-ava-ont-ul"
ASM_OPTS="-h5000 -c2"
MIN_MATCH_FRAC="0.2"
THREADS=4
TIMEOUT="60m"

while getopts "d:Ra:P:r:q:p:m:b:B:H:e:E:J:M:t:v" opt; do
    case "${opt}" in
        d) DIR="${OPTARG}" ;;
        R) RERUN=true ;;
        a) BINNING_BED="${OPTARG}" ;;
        P) BED_GREP_PATTERN="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        q) READS="${OPTARG}" ;;
        p) PAF="${OPTARG}" ;;
        m) MODE="${OPTARG}" ;;
        b) MINIMAP2="${OPTARG}" ;;
        B) MINIASM="${OPTARG}" ;;
        B) HIFIASM="${OPTARG}" ;;
        e) AVA_OPTS="${OPTARG}" ;;
        E) ASM_OPTS="${OPTARG}" ;;
        F) MIN_MATCH_FRAC="${OPTARG}" ;;
        J) JOIN="${OPTARG}" ;;
        M) MARGIN="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        T) TIMEOUT="${TIMEOUT}" ;;
        v) set -x
    esac
done

function log() {
    echo "[$(date +"%Y-%m-%dT%H:%M:%S%z")] $1" >&2
}

function err_exit() {
    log "E: $1"
    exit 1
}

function check_command() {
    if ! command -v $1 >/dev/null; then err_exit "required tool $1 is not found."; fi
}

function run_command() {
    set +e
    timeout "${TIMEOUT}" $@ &
    pid=$!
    set -e

    echo "1000" > "/proc/${pid}/oom_score_adj"
    wait ${pid}
}

# run_command creates a background process
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

# assembler dispatchers
# every dispatcher takes two arguments, input fastq and output gfa filenames. any assembler can be added by wrapping it as a fastq -> gfa converter.
function run_minimap2() {
    if [[ ! -e "${2%.gfa}.ava.done" ]]; then
        log "start minimap2" 2>"${2%.gfa}.log"
        run_command ${MINIMAP2} -t "${THREADS}" ${AVA_OPTS} "$1" "$1"                                                                       2>>"${2%.gfa}.log" > "${2%.gfa}.paf"
        log "finished minimap2" 2>>"${2%.gfa}.log"
        touch "${2%.gfa}.ava.done"
    fi
    if [[ ! -e "${2%.gfa}.asm.done" ]]; then
        log "start miniasm" 2>>"${2%.gfa}.log"
        run_command ${MINIASM} -f "$1" ${ASM_OPTS} <(cat "${2%.gfa}.paf" | awk -vTH="${MIN_MATCH_FRAC}" '$10 / $11 > TH') "${2%.gfa}.paf"   2>>"${2%.gfa}.log" > "$2"
        log "finished miniasm" 2>>"${2%.gfa}.log"
        touch "${2%.gfa}.asm.done"
    fi
}
function run_hifiasm() {
    if [[ ! -e "${2%.gfa}.asm.done" ]]; then
        log "start hifiasm" 2>"${2%.gfa}.log"
        run_command ${HIFIASM} -t "${THREADS}" "${ASM_OPTS}" -o "${2%.gfa}" "$1"                                                            2>>"${2%.gfa}.log"
        ln -srf "${2%.gfa}.bp.p_ctg.gfa" "$2"
        log "finished hifiasm" 2>>"${2%.gfa}.log"
        touch "${2%.gfa}.asm.done"
    fi
}
function run_asm() {
    case ${MODE} in
        minimap2)   run_minimap2 $1 $2 | timeout "${TIMEOUT}" cat || true ;;
        hifiasm)    run_hifiasm $1 $2 | timeout "${TIMEOUT}" cat || true ;;
    esac
}

# make sure input files exist
[[ -z "${BINNING_BED}" ]]                           && err_exit "-a <binning_bed> is required."
[[ ! -f "${BINNING_BED}" ]]                         && err_exit "-a ${BINNING_BED} does not exist or is not a regular file."
[[ ! "${BINNING_BED}" =~ .+\.bed$ ]]                && err_exit "-a ${BINNING_BED} must be in the BED format."
BINNING_BED=$(realpath "${BINNING_BED}")

[[ -z "${REF}" ]]                                   && err_exit "-r <reference> is required."
[[ ! -f "${REF}" ]]                                 && err_exit "-r ${REF} does not exist or is not a regular file."
[[ ! "${REF}" =~ .+\.(fa|fasta|fa.gz|fasta.gz)$ ]]  && err_exit "-r ${REF} must be in the FASTA format."
REF=$(realpath "${REF}")

if [[ -n "${PAF}" ]]; then
    [[ ! -f "${PAF}" ]]                             && err_exit "-p ${PAF} does not exist or is not a regular file."
    [[ ! "${PAF}" =~ .+\.paf$ ]]                    && err_exit "-p ${PAF} must be in the PAF format."
    PAF=$(realpath "${PAF}")
fi

[[ -z "${READS}" ]]                                 && err_exit "-q <reads> is required."
[[ ! -f "${READS}" ]]                               && err_exit "-q ${READS} does not exist or is not a regular file."
[[ ! "${READS}" =~ .+\.(fq|fastq|fq.gz|fastq.gz) ]] && err_exit "-q ${READS} must be in the FASTQ format."
READS=$(realpath "${READS}")

# check mode
case "${MODE}" in
    minimap2|hifiasm) ;;
    *) err_exit "-m <mode> must be one of minimap2 or hifiasm." ;;
esac

# check depending commands exist
check_command "zcat"
check_command "seqkit"
check_command "bedtools"
[[ -z "${PAF}" ]]                                   && check_command "${MINIMAP2}"  # for generating paf
[[ ! -f "${REF}.fai" ]]                             && check_command "samtools"     # for generating index
[[ "${MODE}" = "minimap2" ]]                        && check_command "${MINIMAP2}" && check_command "${MINIASM}"
[[ "${MODE}" = "hifiasm" ]]                         && check_command "${HIFIASM}"

# determine in which directory to run
if ${RERUN}; then
    [[ -z "${DIR}" ]]                               && err_exit "-d <dir> is required when rerun."
    [[ ! -d "${DIR}" ]]                             && err_exit "-d ${DIR} does not exist or is not a directory."
    DIR=$(realpath "${DIR}")
elif [[ -n "${DIR}" ]]; then
    [[ -e "${DIR}" ]]                               && err_exit "-d ${DIR} already exists."
else
    DIR="$(pwd)/run_${MODE}_$(date +%Y%m%d%H%M)"
fi
mkdir -p "${DIR}"
pushd "${DIR}"
echo "$@" > args.txt

# create index of reference if it doesn't exist
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi

# filter annotations
log "filtering annotations"
cat "${BINNING_BED}" | awk -vp="${BED_GREP_PATTERN}" '$4 ~ p' > "annot.flt.bed"
bedtools merge -i "annot.flt.bed" -d "${JOIN}" -c 4,5,6 -o distinct,distinct,distinct > "annot.flt.merged.bed"
bedtools slop -i "annot.flt.merged.bed" -g "${REF}.fai" -b "${MARGIN}" > "annot.flt.merged.margined.bed"
log "finished filtering annotations"

# create dirs for each regions
log "creating directories for each regions"
cat "annot.flt.merged.margined.bed" | awk '{ printf("%s_%s_%s\n", $1, $2, $3); }' | xargs mkdir -p
cat "annot.flt.merged.margined.bed" | awk '{ file=sprintf("%s_%s_%s/region.bed", $1, $2, $3); print $0 > file }'

# create paf if it doesn't exist
if [[ ! -e "${PAF}" ]]; then
    log "mapping ${READS} to ${REF} -> ${PAF}"
    READS_BASE=$(basename "${READS%.(fq|fastq|fq.gz|fastq.gz)}")
    REF_BASE=$(basename "${REF%.(fa|fasta|fa.gz|fasta.gz)}")
    PAF="${READS_BASE}.${REF_BASE}.xmapont.paf"
    minimap2 -t"${THREADS}" -xmap-ont "${REF}" "${READS}" > "${PAF}"
    log "mapping finished"
else
    log "mapping already exists: ${PAF}"
fi

# build read name -> region map
if [[ ! -e "names.sorted.joined.txt" ]]; then
    log "generating read-to-region map"
    cat "${PAF}" | grep "tp:A:P" | awk -vm="${MIN_SPAN}" '$4 - $3 > m' > "filtered.paf"
    for r in $(ls ./*/region.bed); do
        TAG=$(basename $(dirname "$r"))
        CHR=$(cat "$r" | cut -f1)
        BEG=$(cat "$r" | cut -f2)
        END=$(cat "$r" | cut -f3)
        cat "filtered.paf" | awk -vt="${TAG}" -vc="${CHR}" -vb="${BEG}" -ve="${END}" '{ if ($6 == c && $8 < e && $9 >= b) printf "%s\t%s\n", t, $1; }'
    done > "names.txt"

    cat "names.txt" | sort -k2,2 -k1,1 | uniq | uniq -f1 --group=append \
        | awk '{ if ($0 == "") { printf "%s%s\n", n, a; n=""; a=""; } else { n=$2; a=sprintf("%s\t%s", a, $1); } }' > "names.sorted.joined.txt"
else
    log "read-to-region map already exists"
fi

# sort the map to the order the reads are stored in the input
if [[ ! -e "names.sorted.joined.ordered.txt" ]]; then
    log "collecting read name order in ${READS}"
    seqkit fx2tab "${READS}" | awk '{ printf "%s\t%09d\n", $1, NR }' | sort > "order.txt"
    join -11 -21 -t$'\t' "order.txt" "names.sorted.joined.txt" | cut -f2- | sort > "names.sorted.joined.ordered.txt"
    log "finished collecting read name order"
else
    log "read name order in ${READS} already exists"
fi

# then scatter reads to regions
if [[ ! -e "collect.done" ]]; then
    log "collecting reads"
    join -11 -21 -t$'\t' <(seqkit fx2tab "${READS}" | awk '{ printf "%09d\t%s\n", NR, $0; }') "names.sorted.joined.ordered.txt" \
        | awk '{ r=sprintf("@%s\n%s\n+\n%s", $2, $3, $4); for (i = 5; i <= NF; i++) { f=sprintf("%s/reads.fq", $i); print r >> f } }'
    touch "collect.done"
    log "finished collecting reads"
else
    log "reads are already collected"
fi

# all-vs-all -> asm
if ! [[ -e "contigs.fa" ]]; then
    for r in $(cat "annot.flt.merged.margined.bed" | awk '{ printf("%s_%s_%s\n", $1, $2, $3); }'); do
        log "running assembler on ${r}"
        run_asm "${r}/reads.fq" "${r}/asm.gfa"
        cat "${r}/asm.gfa" | (grep "^S" || true) | awk -vr="${r}" '{ printf ">%s_%s\n%s\n", r, $2, $3; }' > "${r}/contigs.fa"
    done
    log "finished assembly"

    for r in $(cat "annot.flt.merged.margined.bed" | awk '{ printf("%s_%s_%s\n", $1, $2, $3); }'); do
        cat "${r}/contigs.fa"
    done > "contigs.fa"
else
    log "contigs.fa already exist"
fi

popd

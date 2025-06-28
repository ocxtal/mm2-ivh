#! /bin/bash
set -eux -o pipefail

SCRIPT_DIR="$(dirname "$(realpath "$0")")"

bash "${SCRIPT_DIR}/install_commands.sh"
bash "${SCRIPT_DIR}/download_data.sh"

samtools view -S -F2816 data/all_pass.vhg002v1.bam | awk '{ printf("@%s\n%s\n+\n%s\n", $1, $10, $11); }' > all.fq
filtlong --min_length 50000 --min_mean_q 90.0 all.fq > all.flt.m50k.90.fq
minimap2 -t32 -xmap-ont data/chm13v2.0.fa.gz all.flt.m50k.90.fq > all.flt.m50k.90.chm13.v2.0.xmapont.paf

# mm2-ivh
bash "${SCRIPT_DIR}/run_local_assembly.sh" \
    -m minimap2 \
    -b "$(realpath bin/mm2-ivh)" \
    -e "-xivh-ava-ont-ul -H -k19 -w15 --wing=3" \
    -E "-h5000 -c2" \
    -d mm2-ivh-hpc \
    -t 32 \
    -r "$(realpath data/chm13v2.0.fa.gz)" \
    -a "$(realpath data/chm13v2.0_censat_v2.1.bed)" \
    -p all.flt.m50k.90.chm13.v2.0.xmapont.paf \
    -q "$(realpath all.flt.m50k.90.fq)"

# mm2-k133-hpc
bash "${SCRIPT_DIR}/run_local_assembly.sh" \
    -m minimap2 \
    -b "$(realpath bin/mm2-ivh)" \
    -e "-xivh-ava-ont-ul -H -k133 -w15 --wing=0" \
    -E "-h5000 -c2" \
    -d mm2-k133-hpc \
    -t 32 \
    -r "$(realpath data/chm13v2.0.fa.gz)" \
    -a "$(realpath data/chm13v2.0_censat_v2.1.bed)" \
    -p all.flt.m50k.90.chm13.v2.0.xmapont.paf \
    -q "$(realpath all.flt.m50k.90.fq)"

# hifiasm-ont-r1
bash "${SCRIPT_DIR}/run_local_assembly.sh" \
    -m hifiasm \
    -b "$(realpath bin/hifiasm)" \
    -e "--ont -r1" \
    -d hifiasm-ont-r1 \
    -t 32 \
    -r "$(realpath data/chm13v2.0.fa.gz)" \
    -a "$(realpath data/chm13v2.0_censat_v2.1.bed)" \
    -p all.flt.m50k.90.chm13.v2.0.xmapont.paf \
    -q "$(realpath all.flt.m50k.90.fq)"

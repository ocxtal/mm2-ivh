#! /bin/bash
set -eu -o pipefail

DIR="data"

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

# check prereqs
check_command "wget"
check_command "md5sum"
check_command "aws"

if [[ ! -d "${DIR}" ]]; then
  echo "Directory ${DIR} does not exist, creating it."
  mkdir -p "${DIR}"
fi

# CHM13 v2.0
CHM13_GENOME_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
CHM13_GENOME_MD5="06caded0055ec647a69acb708e13beff"

if [[ ! -e "${DIR}/chm13v2.0.fa.gz" ]]; then
  echo "Downloading CHM13 v2.0 genome."
  wget "${CHM13_GENOME_URL}" -O "${DIR}/chm13v2.0.fa.gz"
else
  echo "CHM13 v2.0 genome already downloaded."
fi

echo "Checking MD5 checksum for CHM13 v2.0 genome."
echo "${CHM13_GENOME_MD5}  ${DIR}/chm13v2.0.fa.gz" | md5sum -c --status

# CHM13 v2.0 annotation
CHM13_ANNOT_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed"
CHM13_ANNOT_MD5="9838d9f172e34930425d7582243d4041"

if [[ ! -e "${DIR}/chm13v2.0_censat_v2.1.bed" ]]; then
  echo "Downloading CHM13 v2.0 annotation."
  wget "${CHM13_ANNOT_URL}" -O "${DIR}/chm13v2.0_censat_v2.1.bed"
else
  echo "CHM13 v2.0 annotation already downloaded."
fi

echo "Checking MD5 checksum for CHM13 v2.0 annotation."
echo "${CHM13_ANNOT_MD5}  ${DIR}/chm13v2.0_censat_v2.1.bed" | md5sum -c --status

# HG002 v1.1
HG002_GENOME_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz"
HG002_GENOME_MD5="151b2eb3657cc5c3559e7cecb86303f4"

if [[ ! -e "${DIR}/hg002v1.1.fasta.gz" ]]; then
  echo "Downloading HG002 v1.1 genome."
  wget "${HG002_GENOME_URL}" -O "${DIR}/hg002v1.1.fasta.gz"
else
  echo "HG002 v1.1 genome already downloaded."
fi

echo "Checking MD5 checksum for HG002 v1.1 genome."
echo "${HG002_GENOME_MD5}  ${DIR}/hg002v1.1.fasta.gz" | md5sum -c --status

# HG002 v1.1 annotation
HG002_ANNOT_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.bed"
HG002_ANNOT_MD5="c34dbed7002e8a07b53c79b08237ebbb"

if [[ ! -e "${DIR}/hg002v1.1.cenSatv2.0.bed" ]]; then
  echo "Downloading HG002 v1.1 annotation."
  wget "${HG002_ANNOT_URL}" -O "${DIR}/hg002v1.1.cenSatv2.0.bed"
else
  echo "HG002 v1.1 annotation already downloaded."
fi

echo "Checking MD5 checksum for HG002 v1.1 annotation."
echo "${HG002_ANNOT_MD5}  ${DIR}/hg002v1.1.cenSatv2.0.bed" | md5sum -c --status

# HG002 reads
HG002_READS_URL="s3://ont-open-data/gm24385_2023.12/all_pass.vhg002v1.bam"
HG002_READS_MD5="014f635ecbb27d4992ae9a9580775bfa"

if [[ ! -e "${DIR}/all_pass.vhg002v1.bam" ]]; then
  echo "Downloading HG002 reads."
  aws s3 --no-sign-request cp "${HG002_READS_URL}" "${DIR}/all_pass.vhg002v1.bam"
else
  echo "HG002 reads already downloaded."
fi

echo "Checking MD5 checksum for HG002 reads."
echo "${HG002_READS_MD5}  ${DIR}/all_pass.vhg002v1.bam" | md5sum -c --status

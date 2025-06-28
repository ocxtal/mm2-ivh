#! /bin/bash
set -eu -o pipefail

ROOT_DIR="."
SRC_DIR="${ROOT_DIR}/src"
BIN_DIR="${ROOT_DIR}/bin"
REPO_DIR="$(git rev-parse --show-toplevel)"

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

check_command "wget"

if [[ ! -d "${SRC_DIR}" ]]; then
    mkdir -p "${SRC_DIR}"
fi

if [[ ! -d "${BIN_DIR}" ]]; then
    mkdir -p "${BIN_DIR}"
fi

# mm2-ivh
MM2_IVH_VERSION="HEAD"
if [[ ! -e "${BIN_DIR}/mm2-ivh" ]]; then
    echo "Installing mm2-ivh."
    if [[ ! -d "${SRC_DIR}/mm2-ivh" ]]; then
        git clone "${REPO_DIR}" "${SRC_DIR}/mm2-ivh"
    fi
    pushd "${SRC_DIR}/mm2-ivh"
    git fetch origin --tags
    git checkout "${MM2_IVH_VERSION}"
    make -j
    popd
    cp "${SRC_DIR}/mm2-ivh/mm2-ivh" "${BIN_DIR}/mm2-ivh"
else
    echo "mm2-ivh is already installed."
fi

# miniasm
MINIASM_VERSION="ce615d1d6b8678d38f2f9d27c9dccd944436ae75"
if [[ ! -e "${BIN_DIR}/miniasm" ]]; then
    echo "Installing miniasm."
    if [[ ! -d "${SRC_DIR}/miniasm" ]]; then
        git clone https://github.com/lh3/miniasm.git "${SRC_DIR}/miniasm"
    fi
    pushd "${SRC_DIR}/miniasm"
    git fetch origin --tags
    git checkout "${MINIASM_VERSION}"
    make -j
    popd
    cp "${SRC_DIR}/miniasm/miniasm" "${BIN_DIR}/miniasm"
else
    echo "miniasm is already installed."
fi

# hifiasm
HIFIASM_VERSION="0.25.0"
if [[ ! -e "${BIN_DIR}/hifiasm" ]]; then
    echo "Installing hifiasm."
    if [[ ! -d "${SRC_DIR}/hifiasm" ]]; then
        git clone https://github.com/chhylp123/hifiasm.git "${SRC_DIR}/hifiasm"
    fi
    pushd "${SRC_DIR}/hifiasm"
    git fetch origin --tags
    git checkout "${HIFIASM_VERSION}"
    make -j
    popd
    cp "${SRC_DIR}/hifiasm/hifiasm" "${BIN_DIR}/hifiasm"
else
    echo "hifiasm is already installed."
fi

# seqkit
SEQKIT_VERSION="2.8.2"
SEQKIT_OS=$(uname -s | tr '[:upper:]' '[:lower:]')
SEQKIT_ARCH=$(uname -m | sed 's/x86_64/amd64/' | sed 's/aarch64/arm64/')
if [[ ! -e "${BIN_DIR}/seqkit" ]]; then
    echo "Installing seqkit."
    wget "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_${SEQKIT_OS}_${SEQKIT_ARCH}.tar.gz" -O "${SRC_DIR}/seqkit.tar.gz"
    tar -xzf "${SRC_DIR}/seqkit.tar.gz" -C "${BIN_DIR}"
else
    echo "seqkit is already installed."
fi    

# bedtools
BEDTOOLS_VERSION="v2.31.1"
if [[ ! -e "${BIN_DIR}/bedtools" ]]; then
    echo "Installing bedtools."
    if [[ ! -d "${SRC_DIR}/bedtools2" ]]; then
        git clone https://github.com/arq5x/bedtools2.git "${SRC_DIR}/bedtools2"
    fi
    pushd "${SRC_DIR}/bedtools2"
    git checkout "${BEDTOOLS_VERSION}"
    make -j
    popd
    find "${SRC_DIR}/bedtools2/bin" -maxdepth 1 -type f -executable -exec cp {} "${BIN_DIR}/." \;
else
    echo "bedtools is already installed."
fi

# samtools
SAMTOOLS_VERSION="1.20"
if [[ ! -e "${BIN_DIR}/samtools" ]]; then
    echo "Installing samtools."
    if [[ ! -d "${SRC_DIR}/samtools" ]]; then
        git clone https://github.com/samtools/samtools.git "${SRC_DIR}/samtools"
    fi
    if [[ ! -d "${SRC_DIR}/htslib" ]]; then
        git clone https://github.com/samtools/htslib.git "${SRC_DIR}/htslib"
    fi

    pushd "${SRC_DIR}/htslib"
    git fetch origin --tags
    git submodule update --init --recursive
    git checkout "${SAMTOOLS_VERSION}"
    make -j
    popd

    pushd "${SRC_DIR}/samtools"
    git fetch origin --tags
    git submodule update --init --recursive
    git checkout "${SAMTOOLS_VERSION}"
    make -j
    popd

    cp "${SRC_DIR}/samtools/samtools" "${BIN_DIR}/samtools"
else
    echo "samtools is already installed."
fi

# filtlong
FILTLONG_VERSION="7c654f1df394fbe745e755357c10081860bfa5ae"
if [[ ! -e "${BIN_DIR}/filtlong" ]]; then
    echo "Installing filtlong."
    if [[ ! -d "${SRC_DIR}/Filtlong" ]]; then
        git clone https://github.com/rrwick/Filtlong.git "${SRC_DIR}/Filtlong"
    fi
    pushd "${SRC_DIR}/Filtlong"
    git checkout "${FILTLONG_VERSION}"
    make -j
    popd
    cp "${SRC_DIR}/Filtlong/filtlong" "${BIN_DIR}/filtlong"
else
    echo "filtlong is already installed."
fi

# minimap2
MINIMAP2_VERSION="2.28"
if [[ ! -e "${BIN_DIR}/minimap2" ]]; then
    echo "Installing minimap2."
    if [[ ! -d "${SRC_DIR}/minimap2" ]]; then
        git clone https://github.com/lh3/minimap2.git "${SRC_DIR}/minimap2"
    fi
    pushd "${SRC_DIR}/minimap2"
    git fetch origin --tags
    git checkout "${MINIMAP2_VERSION}"
    make -j
    popd
    cp "${SRC_DIR}/minimap2/minimap2" "${BIN_DIR}/minimap2"
else
    echo "minimap2 is already installed."
fi

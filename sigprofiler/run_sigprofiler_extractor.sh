#!/bin/bash
#SBATCH -t 1-0:0:0
#SBATCH -c 1
#SBATCH --mem 160G
#SBATCH -p superhimem
module load python3
# run_sigprofiler_extractor.sh
# Runs SigProfilerExtractor (de novo SBS/DBS/ID signature extraction) over the
# same cohort of VEP-annotated VarScan2 somatic VCFs used by
# parse_vcf_varscan_vep.py:
#   $INPUT_DIR/{folder}/{sample}.wgs.VarScan2.somatic_annotation.vcf.gz
#
# The Python driver decompresses every *.vcf.gz into a single flat staging
# directory (SigProfilerMatrixGenerator requires uncompressed, one-per-sample
# VCFs in one folder), then extracts and decomposes signatures.
#
# One-time setup (uncomment to install, or run --install-genome once):
#   pip install SigProfilerExtractor
#   # the --install-genome flag below downloads GRCh38 on the first run
#
# Usage:
#   bash run_sigprofiler_extractor.sh
#
# Output: SigProfilerExtractor results tree under $OUT_DIR
#   (SBS96/DBS78/ID83 matrices, De_Novo_Solution/, Suggested_Solution/,
#    COSMIC decomposition, and selection plots)

set -euo pipefail

INPUT_DIR="/cluster/projects/vannergroup/controlled_data/GDC/2026-06-30"
SCRIPT="/cluster/projects/vannergroup/Stueckmann/scripts/sigProfiler/run_sigprofiler_extractor.py"
OUT_DIR="/cluster/projects/vannergroup/Stueckmann/out/sigProfiler"


GENOME="GRCh38"
MIN_SIG=1
MAX_SIG=10
REPLICATES=100
CPU=-1   # -1 = use all available cores

# Optional: store/read the ~3 GB reference genome outside the installed
# package (e.g. a shared, writable project dir). Leave empty to use the
# default location inside the SigProfilerMatrixGenerator package.
# NOTE: SigProfiler resolves the genome at <VOLUME>/tsb/<genome>, so the
# GRCh38 chromosome files must live at <VOLUME>/tsb/GRCh38/ (not <VOLUME>/GRCh38).
VOLUME="/cluster/projects/vannergroup/Stueckmann/input/sigProfiler"
VOLUME_FLAG=""
if [ -n "$VOLUME" ]; then
    mkdir -p "$VOLUME"
    VOLUME_FLAG="--volume $VOLUME"
fi

mkdir -p "$OUT_DIR"

# Set INSTALL_GENOME=1 on the first run to download the GRCh38 reference
# (~3 GB, one time per machine); leave 0 afterwards.
INSTALL_GENOME=0
INSTALL_FLAG=""
if [ "$INSTALL_GENOME" = "1" ]; then
    INSTALL_FLAG="--install-genome"
fi

python3 "$SCRIPT" \
    "$INPUT_DIR" \
    "$OUT_DIR" \
    --genome "$GENOME" \
    --min-signatures "$MIN_SIG" \
    --max-signatures "$MAX_SIG" \
    --replicates "$REPLICATES" \
    --cpu "$CPU" \
    --pass-only \
    $VOLUME_FLAG \
    $INSTALL_FLAG

echo ""
echo "Done. Signature results under: $OUT_DIR"

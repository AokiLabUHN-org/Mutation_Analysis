#!/bin/bash
#SBATCH -t 5-0:0:0
#SBATCH -c 1
#SBATCH --mem 30G
# run_parse_vcf_varscan_vep.sh
# Runs parse_vcf_varscan_vep.py on every
# *.wgs.VarScan2.somatic_annotation.vcf.gz found under the GDC
# controlled-data directory (~11,000 files, one per sample subfolder),
# merging all per-sample results into a single TSV with a sample_id column.
#
# Usage:
#   bash run_parse_vcf_varscan_vep.sh
#
# Output: $MERGED_TSV (one row per variant, all samples combined)
# Errors: $ERROR_LOG (files that failed to parse; run continues past them)

set -euo pipefail

INPUT_DIR="/cluster/projects/vannergroup/controlled_data/GDC/2026-06-30"
SCRIPT="/cluster/projects/vannergroup/Stueckmann/scripts/dndscv/parse_vcf_varscan_vep.py"
OUT_DIR="/cluster/projects/vannergroup/Stueckmann/out/dndscv/parsed_vcfs"
SUFFIX=".wgs.VarScan2.somatic_annotation.vcf.gz"

MERGED_TSV="${OUT_DIR}/GDC_VarScan2.merged.parsed.tsv"
ERROR_LOG="${OUT_DIR}/GDC_VarScan2.parse_errors.log"

mkdir -p "$OUT_DIR"
: > "$ERROR_LOG"

TMP_TSV=$(mktemp)
trap 'rm -f "$TMP_TSV"' EXIT

n=0
n_ok=0
n_fail=0
header_written=false

# Files live one-per-subfolder: $INPUT_DIR/{folder}/{sample}.wgs.VarScan2.somatic_annotation.vcf.gz
# Process substitution (not a pipe) keeps the loop in the current shell, so
# the counters below survive across iterations.
while read -r vcf; do

    # Derive sample ID from filename: <sample_id>.wgs.VarScan2.somatic_annotation.vcf.gz -> <sample_id>
    filename=$(basename "$vcf")
    sample="${filename%"$SUFFIX"}"

    n=$((n + 1))
    echo "[$(date '+%H:%M:%S')] (${n}) Processing ${sample}..."

    if python3 "$SCRIPT" "$vcf" "$TMP_TSV" > /dev/null 2>> "$ERROR_LOG"; then
        n_ok=$((n_ok + 1))
    else
        echo "  WARNING: failed to parse ${sample}, see ${ERROR_LOG}"
        n_fail=$((n_fail + 1))
        continue
    fi

    if [ "$header_written" = false ]; then
        head -n 1 "$TMP_TSV" | awk -F'\t' -v OFS='\t' '{print "sample_id", $0}' > "$MERGED_TSV"
        header_written=true
    fi

    # Append data rows with sample_id prepended
    tail -n +2 "$TMP_TSV" | awk -F'\t' -v OFS='\t' -v s="$sample" '{print s, $0}' >> "$MERGED_TSV"

done < <(find "$INPUT_DIR" -mindepth 2 -maxdepth 2 -name "*${SUFFIX}" | sort)

echo ""
echo "Done."
echo "Files found:   ${n}"
echo "Parsed OK:     ${n_ok}"
echo "Failed:        ${n_fail} (see ${ERROR_LOG})"
echo "Merged output: ${MERGED_TSV}"
wc -l "$MERGED_TSV"

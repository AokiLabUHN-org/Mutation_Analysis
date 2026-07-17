#!/usr/bin/env python3
"""
parse_vcf_varscan_vep.py
Parses a VEP-annotated VarScan2 somatic VCF and extracts per-variant fields.

Adapted from parse_vcf.py (which targets Funcotator-annotated Mutect2 VCFs)
to match this pipeline's actual VCF format:
  - Consequence annotation lives in the INFO/CSQ field (VEP), not FUNCOTATION
  - Sample FORMAT is VarScan2-style: GT:GQ:DP:RD:AD:FREQ:DP4
    (RD = ref reads, AD = alt reads, FREQ = allele frequency as a percent)
  - Two sample columns are present (NORMAL, TUMOR); the TUMOR column is used

Output columns:
  gene             VEP SYMBOL
  location         CHROM:POS
  chrom / pos      Split location fields
  ref / alt        Reference and alternate alleles
  mut_type         SNV / insertion / deletion / MNV / indel
  VAF              Allele fraction (FORMAT FREQ, converted to 0-1 fraction)
  alt_reads        Alt-supporting reads (FORMAT AD)
  ref_reads        Ref reads (FORMAT RD)
  total_depth      Total depth (FORMAT DP)
  consequence      Simplified consequence (missense/nonsense/silent/etc.)
  raw_consequence  Full VEP Consequence string
  FILTER           VCF FILTER field

Usage:
  python parse_vcf_varscan_vep.py input.vep.vcf[.gz] [output.tsv]
"""

import sys
import csv
import re


# ---------------------------------------------------------------------------
# Mutation type classification
# ---------------------------------------------------------------------------

def classify_mutation_type(ref, alts):
    types = []
    for alt in alts:
        if len(ref) == 1 and len(alt) == 1:
            types.append("SNV")
        elif len(ref) == 1 and len(alt) > 1:
            types.append("insertion")
        elif len(ref) > 1 and len(alt) == 1:
            types.append("deletion")
        elif len(ref) != len(alt):
            types.append("indel")
        else:
            types.append("MNV")
    return ",".join(types)


# ---------------------------------------------------------------------------
# FORMAT/SAMPLE field parser
# ---------------------------------------------------------------------------

def parse_format(fmt_str, sample_str):
    return dict(zip(fmt_str.split(":"), sample_str.split(":")))


# ---------------------------------------------------------------------------
# VEP CSQ field name extraction from ##INFO header
# ---------------------------------------------------------------------------

def parse_csq_fields(meta_lines):
    """Return ordered list of VEP CSQ field names from the ##INFO header."""
    for line in meta_lines:
        if "##INFO=<ID=CSQ" in line:
            m = re.search(r'Format: ([^"]+)"', line)
            if m:
                return [f.strip() for f in m.group(1).split("|")]
    return None


# ---------------------------------------------------------------------------
# CSQ INFO field parser
#
# The annotation is encoded as: CSQ=allele1_fields,allele2_fields,...
# where each comma-separated block is one transcript annotation with
# pipe-delimited fields. Takes the first transcript annotation if
# multiple are present.
# ---------------------------------------------------------------------------

def extract_csq(info_str, field_names):
    """
    Parse CSQ from INFO string.
    Returns (gene, raw_consequence).
    Takes the first transcript annotation if multiple are present.
    """
    csq_val = None
    for kv in info_str.split(";"):
        if kv.startswith("CSQ="):
            csq_val = kv[len("CSQ="):]
            break

    if not csq_val:
        return "N/A", "N/A"

    # Take first transcript (comma-separated at top level)
    first_tx = csq_val.split(",")[0]
    parts = first_tx.split("|")

    if field_names and len(parts) >= 4:
        fd = dict(zip(field_names, parts))
        gene        = fd.get("SYMBOL", "N/A")
        consequence = fd.get("Consequence", "N/A")
    else:
        # Positional fallback: field 3 = SYMBOL, field 1 = Consequence
        gene        = parts[3] if len(parts) > 3 else "N/A"
        consequence = parts[1] if len(parts) > 1 else "N/A"

    return gene.strip(), consequence.strip()


# ---------------------------------------------------------------------------
# Consequence simplification map (VEP Sequence Ontology terms)
# ---------------------------------------------------------------------------

CONSEQUENCE_MAP = {
    "SYNONYMOUS_VARIANT":              "silent",
    "STOP_RETAINED_VARIANT":           "silent",
    "MISSENSE_VARIANT":                "missense",
    "STOP_GAINED":                     "nonsense",
    "STOP_LOST":                       "stop_lost",
    "START_LOST":                      "start_codon",
    "FRAMESHIFT_VARIANT":              "frameshift",
    "INFRAME_INSERTION":               "inframe_indel",
    "INFRAME_DELETION":                "inframe_indel",
    "SPLICE_DONOR_VARIANT":            "splice_site",
    "SPLICE_ACCEPTOR_VARIANT":         "splice_site",
    "SPLICE_REGION_VARIANT":           "splice_site",
    "SPLICE_DONOR_5TH_BASE_VARIANT":   "splice_site",
    "SPLICE_DONOR_REGION_VARIANT":     "splice_site",
    "SPLICE_POLYPYRIMIDINE_TRACT_VARIANT": "splice_site",
    "MATURE_MIRNA_VARIANT":            "RNA",
    "NON_CODING_TRANSCRIPT_EXON_VARIANT": "RNA",
    "NON_CODING_TRANSCRIPT_VARIANT":   "RNA",
    "INTERGENIC_VARIANT":              "intergenic",
    "INTRON_VARIANT":                  "intronic",
    "5_PRIME_UTR_VARIANT":             "5_prime_UTR",
    "3_PRIME_UTR_VARIANT":             "3_prime_UTR",
    "UPSTREAM_GENE_VARIANT":           "5_prime_flank",
    "DOWNSTREAM_GENE_VARIANT":         "3_prime_flank",
    "PROTEIN_ALTERING_VARIANT":        "protein_altering",
    "CODING_SEQUENCE_VARIANT":         "coding_sequence",
    "INCOMPLETE_TERMINAL_CODON_VARIANT": "start_codon",
}

def simplify_consequence(raw):
    # VEP can report multiple co-occurring terms joined by "&"; classify
    # using the first (most severe, per VEP's own ordering) term.
    first_term = raw.split("&")[0]
    return CONSEQUENCE_MAP.get(first_term.upper(), raw)


# ---------------------------------------------------------------------------
# VAF helper — VarScan2 reports FREQ as a percent string, e.g. "19.4%"
# ---------------------------------------------------------------------------

def freq_to_fraction(freq_str):
    if not freq_str or freq_str == ".":
        return "."
    try:
        return f"{float(freq_str.rstrip('%')) / 100:.4f}"
    except ValueError:
        return freq_str


# ---------------------------------------------------------------------------
# Main VCF parser — expects a single VEP-annotated VarScan2 somatic VCF
# with NORMAL and TUMOR sample columns
# ---------------------------------------------------------------------------

def parse_vcf(vcf_path):
    meta_lines = []
    csq_fields = None
    rows = []
    sample_name = "SAMPLE"
    sample_idx = None

    import gzip
    opener = gzip.open(vcf_path, "rt") if vcf_path.endswith(".gz") else open(vcf_path)

    with opener as fh:
        for line in fh:
            line = line.rstrip("\n")

            # Collect metadata lines
            if line.startswith("##"):
                meta_lines.append(line)
                continue

            # Column header — parse CSQ field names now that we have all
            # ## lines, then locate the TUMOR sample column
            if line.startswith("#CHROM"):
                csq_fields = parse_csq_fields(meta_lines)
                if csq_fields is None:
                    print("ERROR: No CSQ INFO header found. "
                          "Is this a VEP-annotated VCF?", file=sys.stderr)
                    sys.exit(1)
                cols = line.lstrip("#").split("\t")
                for i, c in enumerate(cols):
                    if c.upper() == "TUMOR":
                        sample_idx = i
                        break
                if sample_idx is None and len(cols) > 9:
                    sample_idx = len(cols) - 1
                sample_name = cols[sample_idx] if sample_idx is not None else "SAMPLE"
                continue

            # Data lines
            fields = line.split("\t")
            if len(fields) < 9:
                continue

            chrom   = fields[0]
            pos     = fields[1]
            ref     = fields[3]
            alt_raw = fields[4]
            filt    = fields[6]
            info    = fields[7]
            fmt     = fields[8]
            sample  = fields[sample_idx] if sample_idx is not None and len(fields) > sample_idx else ""

            alts     = alt_raw.split(",")
            mut_type = classify_mutation_type(ref, alts)

            fmt_dict  = parse_format(fmt, sample)
            vaf       = freq_to_fraction(fmt_dict.get("FREQ", "."))
            depth     = fmt_dict.get("DP", ".")
            ref_reads = fmt_dict.get("RD", ".")
            alt_reads = fmt_dict.get("AD", ".")

            gene, raw_csq = extract_csq(info, csq_fields)
            consequence   = simplify_consequence(raw_csq)

            rows.append({
                "gene":            gene,
                "location":        f"{chrom}:{pos}",
                "chrom":           chrom,
                "pos":             pos,
                "ref":             ref,
                "alt":             alt_raw,
                "mut_type":        mut_type,
                "VAF":             vaf,
                "alt_reads":       alt_reads,
                "ref_reads":       ref_reads,
                "total_depth":     depth,
                "consequence":     consequence,
                "raw_consequence": raw_csq,
                "FILTER":          filt,
            })

    return rows, sample_name


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    vcf_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else "parsed_variants.tsv"

    rows, sample_name = parse_vcf(vcf_path)

    fieldnames = [
        "gene",
        "location", "chrom", "pos", "ref", "alt",
        "mut_type",
        "VAF", "alt_reads", "ref_reads", "total_depth",
        "consequence", "raw_consequence",
        "FILTER",
    ]

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Sample:          {sample_name}")
    print(f"Variants parsed: {len(rows)}")
    print(f"Output:          {out_path}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
run_sigprofiler_extractor.py
De novo mutational-signature extraction (SBS / DBS / ID) with
SigProfilerExtractor, run over the SAME cohort of VEP-annotated VarScan2
somatic VCFs consumed by parse_vcf_varscan_vep.py:

    $INPUT_DIR/{folder}/{sample}.wgs.VarScan2.somatic_annotation.vcf.gz
    (~11,000 files, one per sample subfolder)

Why a staging step?
  SigProfilerExtractor's `input_type="vcf"` mode hands the input directory to
  SigProfilerMatrixGenerator, which expects:
    - one *uncompressed* .vcf per sample,
    - all in a single flat directory,
    - the sample name taken from the VCF *filename* (the NORMAL/TUMOR
      genotype columns are ignored — only CHROM/POS/REF/ALT are read).
  Our inputs are gzipped and scattered one-per-subfolder, so we first
  decompress each into a flat staging dir named "<sample_id>.vcf".

Pipeline:
  1. stage_vcfs()   decompress every *.vcf.gz into <staging_dir>/<sample>.vcf
                    (optionally keeping only FILTER==PASS records)
  2. sigProfilerExtractor()  build SBS-96 / DBS-78 / ID-83 matrices, extract
                    de novo signatures across the requested rank range, and
                    decompose them against COSMIC reference signatures.

Prerequisites (one-time):
  pip install SigProfilerExtractor
  # install the reference genome used by the cohort (GDC = GRCh38):
  python -c "from SigProfilerMatrixGenerator import install as g; g.install('GRCh38')"
  # ...or pass --install-genome on the first run below.

Reference genome location:
  By default the genome is stored *inside* the installed
  SigProfilerMatrixGenerator package (references/chromosomes/tsb/<genome>/),
  which is awkward when that path is read-only or shared across envs.
  Two ways to keep it elsewhere:
    (a) --volume /path/to/refgenomes  -> installs to / reads from that dir.
        Requires a SigProfiler version whose install()/sigProfilerExtractor()
        accept a `volume` argument; the script auto-detects support and warns
        (rather than crashing) if the installed version is too old.
    (b) Symlink fallback (works on any version): install once to a shared
        dir, then link it into the package's references directory, e.g.
          pkg=$(python -c 'import SigProfilerMatrixGenerator as m,os;\
                           print(os.path.dirname(m.__file__))')
          ln -s /shared/refgenomes/GRCh38 \
                "$pkg/references/chromosomes/tsb/GRCh38"

Usage:
  python run_sigprofiler_extractor.py INPUT_DIR OUTPUT_DIR [options]

Example:
  python run_sigprofiler_extractor.py \
      /cluster/projects/vannergroup/controlled_data/GDC/2026-06-30 \
      /cluster/projects/vannergroup/Stueckmann/out/dndscv/sigprofiler \
      --pass-only --cpu 16
"""

import argparse
import glob
import gzip
import os
import shutil
import sys

SUFFIX = ".wgs.VarScan2.somatic_annotation.vcf.gz"


# ---------------------------------------------------------------------------
# Staging: gzipped, per-subfolder VCFs -> flat dir of uncompressed <sample>.vcf
# ---------------------------------------------------------------------------

def stage_vcfs(input_dir, staging_dir, suffix=SUFFIX, pass_only=False):
    """
    Decompress every *${suffix} found under input_dir into staging_dir as
    <sample_id>.vcf, where <sample_id> is the filename with the suffix
    stripped (matching parse_vcf_varscan_vep.py's sample-naming).

    Returns the list of sample_ids staged.
    """
    os.makedirs(staging_dir, exist_ok=True)

    # Match the runner's discovery: one file per subfolder, two levels deep.
    pattern = os.path.join(input_dir, "*", "*" + suffix)
    vcf_gzs = sorted(glob.glob(pattern))
    if not vcf_gzs:
        # Fall back to a fully-recursive search in case the layout is deeper.
        pattern = os.path.join(input_dir, "**", "*" + suffix)
        vcf_gzs = sorted(glob.glob(pattern, recursive=True))

    if not vcf_gzs:
        print(f"ERROR: no files matching *{suffix} under {input_dir}",
              file=sys.stderr)
        sys.exit(1)

    samples = []
    for i, vcf_gz in enumerate(vcf_gzs, 1):
        sample = os.path.basename(vcf_gz)[: -len(suffix)]
        out_vcf = os.path.join(staging_dir, sample + ".vcf")
        samples.append(sample)

        # Skip work already done (lets the run resume after an interruption).
        if os.path.exists(out_vcf) and os.path.getsize(out_vcf) > 0:
            continue

        print(f"[stage] ({i}/{len(vcf_gzs)}) {sample}")
        with gzip.open(vcf_gz, "rt") as fin, open(out_vcf, "w") as fout:
            if pass_only:
                # Keep header lines plus only PASS / "." data records, so
                # low-confidence calls don't inflate the signature matrices.
                for line in fin:
                    if line.startswith("#"):
                        fout.write(line)
                        continue
                    cols = line.split("\t")
                    if len(cols) > 6 and cols[6] in ("PASS", "."):
                        fout.write(line)
            else:
                shutil.copyfileobj(fin, fout)

    print(f"[stage] {len(samples)} VCFs ready in {staging_dir}")
    return samples


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def run_extractor(staging_dir, output_dir, genome, min_sig, max_sig,
                  replicates, cpu, context_type, install_genome, volume=None):
    # Imported lazily so --help / staging work without the (heavy) package.
    import inspect
    from SigProfilerExtractor import sigpro as sig

    # By default the reference genome lives inside the installed
    # SigProfilerMatrixGenerator package (references/chromosomes/...). A
    # `volume` points install + extraction at a custom directory instead --
    # useful when the package dir is read-only or shared across envs. The
    # parameter only exists in newer releases, so pass it *only* when the
    # installed function accepts it; otherwise warn rather than crash.
    def supports_volume(fn):
        return "volume" in inspect.signature(fn).parameters

    if install_genome:
        from SigProfilerMatrixGenerator import install as genInstall
        install_kwargs = {}
        if volume:
            if supports_volume(genInstall.install):
                install_kwargs["volume"] = volume
                print(f"[install] genome dir (volume): {volume}")
            else:
                print("[install] WARNING: installed SigProfilerMatrixGenerator "
                      "does not support --volume; installing into the package "
                      "dir. See the symlink fallback in the header docstring.",
                      file=sys.stderr)
        print(f"[install] installing reference genome {genome} "
              "(one-time, ~3 GB download)...")
        genInstall.install(genome, **install_kwargs)

    extract_kwargs = dict(
        input_type="vcf",
        output=output_dir,
        input_data=staging_dir,
        reference_genome=genome,     # genome for the SBS/DBS/ID context calls
        opportunity_genome=genome,   # genome for trinucleotide normalisation
        context_type=context_type,   # "96,DINUC,ID" -> SBS96 + DBS78 + ID83
        minimum_signatures=min_sig,
        maximum_signatures=max_sig,
        nmf_replicates=replicates,
        cpu=cpu,
    )
    if volume:
        if supports_volume(sig.sigProfilerExtractor):
            extract_kwargs["volume"] = volume
            print(f"[extract] genome dir (volume): {volume}")
        else:
            print("[extract] WARNING: installed SigProfilerExtractor does not "
                  "support --volume; expecting the genome inside the package "
                  "dir (or symlinked there). See the header docstring.",
                  file=sys.stderr)

    print(f"[extract] running SigProfilerExtractor on {staging_dir}")
    sig.sigProfilerExtractor(**extract_kwargs)
    print(f"[extract] done. Results under {output_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Mutational-signature extraction with SigProfilerExtractor.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input_dir",
                   help="Root dir containing per-sample subfolders with "
                        f"*{SUFFIX} files.")
    p.add_argument("output_dir",
                   help="Directory for SigProfilerExtractor results.")
    p.add_argument("--staging-dir", default=None,
                   help="Where to write decompressed VCFs "
                        "(default: <output_dir>/staged_vcfs).")
    p.add_argument("--genome", default="GRCh38",
                   help="Reference genome build (GDC data is GRCh38).")
    p.add_argument("--volume", default=None,
                   help="Directory holding the SigProfiler reference genomes, "
                        "instead of the default location inside the installed "
                        "package. Used for both --install-genome and "
                        "extraction (requires a SigProfiler version that "
                        "supports the 'volume' parameter).")
    p.add_argument("--min-signatures", type=int, default=1,
                   help="Minimum number of signatures to test.")
    p.add_argument("--max-signatures", type=int, default=10,
                   help="Maximum number of signatures to test.")
    p.add_argument("--replicates", type=int, default=100,
                   help="NMF replicates per rank.")
    p.add_argument("--cpu", type=int, default=-1,
                   help="CPU cores (-1 = all available).")
    p.add_argument("--context-type", default="96,DINUC,ID",
                   help="Mutation contexts to extract "
                        "(SBS96, DBS78, ID83).")
    p.add_argument("--pass-only", action="store_true",
                   help="Stage only FILTER==PASS (or '.') records.")
    p.add_argument("--install-genome", action="store_true",
                   help="Install the reference genome before extracting "
                        "(needed once per machine).")
    p.add_argument("--stage-only", action="store_true",
                   help="Decompress/stage the VCFs and stop (no extraction).")
    return p.parse_args()


def main():
    args = parse_args()

    staging_dir = args.staging_dir or os.path.join(args.output_dir,
                                                    "staged_vcfs")
    os.makedirs(args.output_dir, exist_ok=True)

    stage_vcfs(args.input_dir, staging_dir, pass_only=args.pass_only)

    if args.stage_only:
        print("[stage-only] skipping extraction.")
        return

    run_extractor(
        staging_dir=staging_dir,
        output_dir=args.output_dir,
        genome=args.genome,
        min_sig=args.min_signatures,
        max_sig=args.max_signatures,
        replicates=args.replicates,
        cpu=args.cpu,
        context_type=args.context_type,
        install_genome=args.install_genome,
        volume=args.volume,
    )


# SigProfilerExtractor uses multiprocessing; guard the entry point so child
# processes re-import this module cleanly rather than re-running main().
if __name__ == "__main__":
    main()

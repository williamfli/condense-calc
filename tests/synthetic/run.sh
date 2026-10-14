#!/bin/bash
# End-to-end test on simulated reads (~1 min on 2 cores). Run from the repo root, in the condseq env:
#   bash tests/synthetic/run.sh OUTDIR [condseq run options, e.g. --qc verbose]
# Then: python tests/synthetic/check_truth.py OUTDIR
set -euo pipefail
OUT=$(realpath -m "$1"); shift
mkdir -p "$OUT"
[ -f "$OUT/input/samples.tsv" ] || python tests/synthetic/simulate.py "$OUT/input" --seed 0
cat > "$OUT/input/config.yaml" <<CFG
samples: samples.tsv
titration_dir: titration
genomes:
  sim: {fasta: ref/sim.fa, bowtie2_index: ref/sim}
cells:
  SIM: {genome: sim, chroms: [chr1, chr2, chrX]}
threads: 2
CFG
python -m condseq run "$OUT/input/config.yaml" -o "$OUT/run" --cores "${SLURM_CPUS_PER_TASK:-2}" "$@"

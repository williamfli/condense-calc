#!/bin/bash
# End-to-end test on simulated reads (~1 min on 2 cores). Run from the repo root, in the condense-calc env:
#   bash tests/synthetic/run.sh OUTDIR [condense-calc run options, e.g. --qc verbose, or --revised]
# Then: python tests/synthetic/check_truth.py OUTDIR. tests/equivalence/run.sh runs it in the default (Park) mode.
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
bowtie2: bowtie2    # this env's, also in Park mode (whose default would install 2.4.1)
bdata: true         # sliding-window areas, compared with Park's bindata_ver3.py by tests/equivalence
CFG
python -m condense_calc run "$OUT/input/config.yaml" -o "$OUT/run" --cores "${SLURM_CPUS_PER_TASK:-2}" "$@"

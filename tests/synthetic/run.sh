#!/bin/bash
# End-to-end test on simulated reads (~1 min on 2 cores). Run from the repo root, in the condense-calc env:
#   bash tests/synthetic/run.sh OUTDIR [condense-calc run options, e.g. --qc verbose, or --revised]
# Then: python tests/synthetic/check_truth.py OUTDIR. tests/equivalence/run.sh runs it in the default (Park) mode.
set -euo pipefail
mkdir -p "$1"; OUT=$(cd "$1" && pwd -P); shift
mkdir -p "$OUT"
[ -f "$OUT/input/samples.tsv" ] || python tests/synthetic/simulate.py "$OUT/input" --seed 0
# chromatin states for the QC panels after ED Fig. 2g-i: four made-up states in 7 kb blocks
python - "$OUT/input/states.bed" <<'PY'
import sys
with open(sys.argv[1], 'w') as f:
    for chrom, size in (('chr1', 300_000), ('chr2', 200_000), ('chrX', 100_000)):
        for k, st in enumerate(range(0, size, 7000)):
            f.write(f"{chrom}\t{st}\t{min(st + 7000, size)}\t{('1_Active', '2_Transcribed', '3_Repressed', '4_Quiescent')[k % 4]}\n")
PY
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
states: states.bed  # QC panels after ED Fig. 2g-i
CFG
python -m condense_calc run "$OUT/input/config.yaml" -o "$OUT/run" --cores "${SLURM_CPUS_PER_TASK:-2}" "$@"

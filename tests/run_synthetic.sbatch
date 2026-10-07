#!/bin/bash
#SBATCH -J condseq_synth
#SBATCH -c 8
#SBATCH -t 0-02:00
#SBATCH -p test
#SBATCH --mem=24G
#SBATCH -o ../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../lab/logs/condseq_pipeline/%x_%j.err
# End-to-end test on simulated reads. Submit from condense-seq-pipeline/, outside an env.
# Output: ../lab/outputs/condense-seq-pipeline/synthetic/{input,run}
set -euo pipefail
module load python/3.12.11-fasrc02
mamba activate condseq
OUT=$(realpath ../lab/outputs)/condense-seq-pipeline/synthetic
mkdir -p $OUT
[ -f $OUT/input/samples.tsv ] || python tests/simulate.py $OUT/input --seed 0
cat > $OUT/input/config.yaml <<CFG
samples: samples.tsv
titration_dir: titration
cells:
  SIM: {genome: sim, chroms: [chr1, chr2, chrX]}
genomes:
  sim: {fasta: ref/sim.fa, bowtie2_index: ref/sim}
trim_to: 50
bin_sizes: [1000, 10000]
align_threads: 4
mscale: {sp: 1}
fit: {bin_size: 10000, min_tnum: 5, model: sigmoid, method: evolution, min_rsq: 0.5}
reads: {mm: 10, min_mapq: 0, dedup: false}
score: {missing_as_nan: false}
CFG
snakemake -s workflow/Snakefile --configfile $OUT/input/config.yaml --directory $OUT/run --cores 8 -p
echo "=== outputs"; find $OUT/run/tables -type f | sort

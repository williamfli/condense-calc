#!/bin/bash
#SBATCH -J condseq_s0_seed1
#SBATCH -c 8
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=9G
#SBATCH --array=0-1
#SBATCH -o ../../../lab/logs/condseq_pipeline/%x_%A_%a.out
#SBATCH -e ../../../lab/logs/condseq_pipeline/%x_%A_%a.err
# Step 0, mechanism test: realign the same reads with bowtie2 --seed 1 (the pipeline used the
# default, 0). The seed only changes which of several equally good hits a multi-mapping read gets.
# Same command as the Snakefile's align rule otherwise. Sizes from the align rule (seff 2026-10-12).
# Submit from condense-seq-pipeline/geo_test/step0/; then score_10kb.sbatch .../geo_step0/seed1.
set -euo pipefail
module load python/3.12.11-fasrc02
mamba activate condseq
RAW=/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/raw/GSE252941
OUT=/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/pipeline_runs/geo_step0/seed1
IDX=$(realpath ../../../lab/data/ref/hg38)/hg38
mkdir -p $OUT/bam $OUT/logs
NAMES=(H1_NCP_PEG_0_2rep H1_NCP_PEG_6_2rep)
RUNS=("SRR31626988 SRR31626989" "SRR31626986 SRR31626987")      # from samples_step0.tsv
name=${NAMES[$SLURM_ARRAY_TASK_ID]}; r1=""; r2=""
for r in ${RUNS[$SLURM_ARRAY_TASK_ID]}; do r1+=${r1:+,}$RAW/${r}_1.fastq.gz; r2+=${r2:+,}$RAW/${r}_2.fastq.gz; done
bowtie2 -p $SLURM_CPUS_PER_TASK --seed 1 --trim-to 50 -x $IDX -1 $r1 -2 $r2 2> $OUT/logs/$name.log \
  | samtools sort -@ 2 -m 1G -o $OUT/bam/$name.bam - && samtools index $OUT/bam/$name.bam

#!/bin/bash
#SBATCH -J condseq_s0_4dnalign
#SBATCH -c 8
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=9G
#SBATCH --array=0-1
#SBATCH -o ../../../lab/logs/condseq_pipeline/%x_%A_%a.out
#SBATCH -e ../../../lab/logs/condseq_pipeline/%x_%A_%a.err
# Step 0, mechanism test 3: realign the same reads with bowtie2 2.5.5 (the pipeline's) against the
# 4DN index from 06a, instead of NCBI's prebuilt index. Tasks: 0 = input, 1 = ID 6.
# Sizes from the align rule (seff 2026-10-12). Submit from condense-seq-pipeline/geo_test/step0/
# after 06a; then score_10kb.sbatch .../geo_step0/idx_4dn.
set -euo pipefail
module load python/3.12.11-fasrc02
mamba activate condseq
NS=/n/netscratch/zhuang_lab/Lab/wfli/condense-seq
NAMES=(H1_NCP_PEG_0_2rep H1_NCP_PEG_6_2rep)
RUNS=("SRR31626988 SRR31626989" "SRR31626986 SRR31626987")
i=$SLURM_ARRAY_TASK_ID; name=${NAMES[$i]}
IDX=$(find $NS/ref_4dn/hg38 -name '*.rev.1.bt2' | head -1); IDX=${IDX%.rev.1.bt2}
OUT=$NS/pipeline_runs/geo_step0/idx_4dn
mkdir -p $OUT/bam $OUT/logs
RAW=$NS/raw/GSE252941
r1=""; r2=""
for r in ${RUNS[$i]}; do r1+=${r1:+,}$RAW/${r}_1.fastq.gz; r2+=${r2:+,}$RAW/${r}_2.fastq.gz; done
echo "$IDX" > $OUT/logs/$name.index
bowtie2 -p $SLURM_CPUS_PER_TASK --trim-to 50 -x $IDX -1 $r1 -2 $r2 2> $OUT/logs/$name.log \
  | samtools sort -@ 2 -m 1G -o $OUT/bam/$name.bam - && samtools index $OUT/bam/$name.bam

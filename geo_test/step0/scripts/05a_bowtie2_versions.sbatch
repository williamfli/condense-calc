#!/bin/bash
#SBATCH -J condseq_s0_bt2ver
#SBATCH -c 8
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=9G
#SBATCH --array=0-3
#SBATCH -o ../../../lab/logs/condseq_pipeline/%x_%A_%a.out
#SBATCH -e ../../../lab/logs/condseq_pipeline/%x_%A_%a.err
# Step 0, mechanism test 2: realign the same reads with older bowtie2 versions (Park's version is
# not recorded). Tasks: 0,1 = bowtie2 2.3.5.1 (input, ID 6); 2,3 = bowtie2 2.4.5 (input, ID 6).
# Same index and options as the pipeline's align rule; sizes from that rule (seff 2026-10-12).
# The bowtie2 envs are throwaway, on netscratch; the first task that needs one creates it.
# Submit from condense-seq-pipeline/geo_test/step0/; then score_10kb.sbatch .../geo_step0/bt2_<version>.
set -euo pipefail
module load python/3.12.11-fasrc02
NS=/n/netscratch/zhuang_lab/Lab/wfli/condense-seq
VERS=(2.3.5.1 2.3.5.1 2.4.5 2.4.5)
NAMES=(H1_NCP_PEG_0_2rep H1_NCP_PEG_6_2rep H1_NCP_PEG_0_2rep H1_NCP_PEG_6_2rep)
RUNS=("SRR31626988 SRR31626989" "SRR31626986 SRR31626987" "SRR31626988 SRR31626989" "SRR31626986 SRR31626987")
i=$SLURM_ARRAY_TASK_ID; v=${VERS[$i]}; name=${NAMES[$i]}
ENV=$NS/envs/bowtie2-$v
mkdir -p $NS/envs
( flock 9; [ -x $ENV/bin/bowtie2 ] || mamba create -y -q -p $ENV -c conda-forge -c bioconda bowtie2=$v samtools=1.21 ) 9>$NS/envs/.lock-$v
mamba activate $ENV
OUT=$NS/pipeline_runs/geo_step0/bt2_$v
mkdir -p $OUT/bam $OUT/logs
RAW=$NS/raw/GSE252941
IDX=$(realpath ../../../lab/data/ref/hg38)/hg38
r1=""; r2=""
for r in ${RUNS[$i]}; do r1+=${r1:+,}$RAW/${r}_1.fastq.gz; r2+=${r2:+,}$RAW/${r}_2.fastq.gz; done
bowtie2 --version | head -1 > $OUT/logs/$name.version
bowtie2 -p $SLURM_CPUS_PER_TASK --trim-to 50 -x $IDX -1 $r1 -2 $r2 2> $OUT/logs/$name.log \
  | samtools sort -@ 2 -m 1G -o $OUT/bam/$name.bam - && samtools index $OUT/bam/$name.bam

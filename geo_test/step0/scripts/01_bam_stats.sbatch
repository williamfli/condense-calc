#!/bin/bash
#SBATCH -J condseq_s0_bamstats
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p shared
#SBATCH --mem=4G
#SBATCH -o ../../../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../../../lab/logs/condseq_pipeline/%x_%j.err
# Step 0, read-level statistics (01_bam_stats.py). Submit from condense-seq-pipeline/geo_test/step0/.
set -euo pipefail
module load python/3.12.11-fasrc02
mamba activate condseq
python scripts/01_bam_stats.py /n/netscratch/zhuang_lab/Lab/wfli/condense-seq/pipeline_runs/geo_step0/bam \
  results ../../../lab/outputs/condense-seq-pipeline/geo_step0

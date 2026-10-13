#!/bin/bash
#SBATCH -J condseq_env
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p test
#SBATCH --mem=4G
# 1 core, 1.0-2.4 GB peak (seff 2026-10-07/08: condseq, condrepro)
#SBATCH -o ../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../lab/logs/condseq_pipeline/%x_%j.err
# Submit from condense-seq-pipeline/, outside an activated env.
# Explicit prefix: `mamba env create -f` ignores envs_dirs in ~/.condarc and installs into home.
module load python/3.12.11-fasrc02
P=/n/holylabs/zhuang_lab/Lab/wfli/conda/envs/condseq
mamba env create -y -p $P -f env/environment.yml
mamba env export -p $P --no-builds | grep -v "^prefix:" | sed "1s/^name: *$/name: condseq/" > env/environment.lock.yml

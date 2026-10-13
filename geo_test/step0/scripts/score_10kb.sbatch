#!/bin/bash
#SBATCH -J condseq_s0_score
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p shared
#SBATCH --mem=4G
#SBATCH -o ../../../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../../../lab/logs/condseq_pipeline/%x_%j.err
# Step 0 realignment tests: 10 kb bincount -> score -> table for the two BAMs in DIR/bam, with the
# same arguments as the Snakefile (bincount 1.4 GB, score < 0.1 GB at 10 kb; seff 2026-10-12).
# Submit from condense-seq-pipeline/geo_test/step0/ after the alignment job:
#   sbatch scripts/score_10kb.sbatch DIR      (DIR = .../geo_step0/seed1, .../geo_step0/bt2_2.4.5, ...)
set -euo pipefail
module load python/3.12.11-fasrc02
mamba activate condseq
export PYTHONPATH=$(realpath ../../..)/condense-seq-pipeline:${PYTHONPATH:-}
OUT=$1
FA=$(realpath ../../../lab/data/ref/hg38)/hg38.fa
TIT=$(realpath ../../../lab/data/park2025/titration)/H1_NCP_PEG.tsv
CHR="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
mkdir -p $OUT/bin $OUT/score $OUT/tables
for n in H1_NCP_PEG_0_2rep H1_NCP_PEG_6_2rep; do
  python -m condseq bincount --bam $OUT/bam/$n.bam --fasta $FA --chroms $CHR --win 10000 --mm 10 -o $OUT/bin/${n}_10kb_bin.gtab.gz
done
python -m condseq score --test $OUT/bin/H1_NCP_PEG_6_2rep_10kb_bin.gtab.gz --control $OUT/bin/H1_NCP_PEG_0_2rep_10kb_bin.gtab.gz \
  --titration $TIT --tnum 6 --chroms $CHR -o $OUT/score/H1_NCP_PEG_6_2rep_10kb_score.gtab.gz
python -m condseq concat --inputs $OUT/score/H1_NCP_PEG_6_2rep_10kb_score.gtab.gz -o $OUT/tables/H1_NCP_PEG_2rep_10kb_score.gtab.gz

#!/bin/bash
#SBATCH -J condseq_ref
#SBATCH -c 16
#SBATCH -t 0-08:00
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -o ../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../lab/logs/condseq_pipeline/%x_%j.err
# Reference genomes + bowtie2 indexes into ../lab/data/ref/ (~4 GB per genome). NOT RUN YET:
# lab storage is ~99% full, so run it when raw reads are coming. Submit from condense-seq-pipeline/.
#
# Park used the 4DN references "4D_hg38" and "4D_mm10" (batch_submit_ver7.py: org_ref). 4DN's
# hg38 is GRCh38 no-alt analysis set; we use NCBI's copy (UCSC chromosome names, prebuilt
# bowtie2 index). For mm10 we use UCSC mm10 primary chromosomes and build the index.
# Check the 4DN portal before relying on mm10: the exact 4DN file wasn't recorded by Park.
set -euo pipefail
module load python/3.12.11-fasrc02
mamba activate condseq
REF=$(realpath ../lab/data)/ref; mkdir -p $REF/hg38 $REF/mm10
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
cd $REF/hg38
wget -q -c $NCBI/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget -q -c $NCBI/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
gunzip -c GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > hg38.fa && samtools faidx hg38.fa
tar -xzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
for f in GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.*.bt2; do
  ln -sf $f hg38.${f#GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.}; done
cd $REF/mm10
wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip -c mm10.fa.gz > mm10.full.fa && samtools faidx mm10.full.fa
samtools faidx mm10.full.fa $(cut -f1 mm10.full.fa.fai | grep -E '^chr([0-9]+|X|Y|M)$') > mm10.fa
samtools faidx mm10.fa && rm mm10.full.fa*
bowtie2-build --threads 16 mm10.fa mm10
ls -la $REF/hg38 $REF/mm10

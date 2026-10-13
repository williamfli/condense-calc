#!/bin/bash
#SBATCH -J condseq_ref
#SBATCH -c 1
#SBATCH -t 0-08:00
#SBATCH -p shared
#SBATCH --mem=8G
#SBATCH -o ../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../lab/logs/condseq_pipeline/%x_%j.err
# hg38: single-threaded, 28 MB peak, 26 min (seff 2026-10-12). mm10 (not measured) runs
# bowtie2-build, which uses all allocated cores: `sbatch -c 8 slurm/10_reference.sbatch mm10` halves
# its few hours of runtime; 8G is the usual estimate for a mammalian index build.
# Reference genomes + bowtie2 indexes into ../lab/data/ref/ (~4-8 GB per genome; archives are
# deleted after extraction). Submit from condense-seq-pipeline/:
#   sbatch slurm/10_reference.sbatch hg38        # or mm10, or both (default)
#
# Park used the 4DN references "4D_hg38" and "4D_mm10" (batch_submit_ver7.py: org_ref). 4DN's
# hg38 is GRCh38 no-alt analysis set; we use NCBI's copy (UCSC chromosome names, prebuilt
# bowtie2 index). For mm10 we use UCSC mm10 primary chromosomes and build the index.
# Check the 4DN portal before relying on mm10: the exact 4DN file wasn't recorded by Park.
set -euo pipefail
module load python/3.12.11-fasrc02
mamba activate condseq
REF=$(realpath ../lab/data)/ref
GENOMES=" ${@:-hg38 mm10} "
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
if [[ "$GENOMES" == *" hg38 "* ]]; then
mkdir -p $REF/hg38 && cd $REF/hg38
wget -q -c $NCBI/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget -q -c $NCBI/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
gunzip -c GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > hg38.fa && samtools faidx hg38.fa
tar -xzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
for f in $(find . -name '*.bt2' ! -name 'hg38.*'); do      # name the index files hg38.{1..4,rev.1,rev.2}.bt2
  ln -sf $f hg38.$(basename $f | grep -oE '(rev\.)?[0-9]\.bt2$'); done
[ $(ls hg38.*.bt2 | wc -l) -eq 6 ] || { echo 'expected 6 hg38 index files'; exit 1; }
rm GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
fi
if [[ "$GENOMES" == *" mm10 "* ]]; then
mkdir -p $REF/mm10 && cd $REF/mm10
wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip -c mm10.fa.gz > mm10.full.fa && samtools faidx mm10.full.fa
samtools faidx mm10.full.fa $(cut -f1 mm10.full.fa.fai | grep -E '^chr([0-9]+|X|Y|M)$') > mm10.fa
samtools faidx mm10.fa && rm mm10.full.fa*
bowtie2-build --threads ${SLURM_CPUS_PER_TASK:-1} mm10.fa mm10 && rm mm10.fa.gz
fi
ls -laL $REF/*/ && du -sh $REF

#!/bin/bash
# Download a genome FASTA (+ .fai) and its bowtie2 index into OUTDIR/<genome>/.
#   get_reference.sh OUTDIR THREADS hg38 [mm10]      (usually via: condseq reference hg38 OUTDIR)
# hg38: NCBI GRCh38 no-alt analysis set with its prebuilt bowtie2 index (~7 GB, ~30 min).
# mm10: UCSC mm10 chr1-19, X, Y, M; the index is built here (a few hours; use more threads).
# Park et al. used the 4DN versions of these references; multi-mapping reads can land differently.
set -euo pipefail
OUT=$(realpath "$1"); THREADS=$2; shift 2
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
for G in "$@"; do
  mkdir -p "$OUT/$G" && cd "$OUT/$G"
  case $G in
  hg38)
    F=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    wget -q -c $NCBI/$F.gz $NCBI/$F.bowtie_index.tar.gz
    gunzip -c $F.gz > hg38.fa && samtools faidx hg38.fa
    tar -xzf $F.bowtie_index.tar.gz
    for f in $(find . -name '*.bt2' ! -name 'hg38.*'); do      # name the index files hg38.*.bt2
      ln -sf "$f" hg38.$(basename "$f" | grep -oE '(rev\.)?[0-9]\.bt2$'); done
    [ "$(ls hg38.*.bt2 | wc -l)" -eq 6 ] || { echo 'expected 6 hg38 index files' >&2; exit 1; }
    rm $F.gz $F.bowtie_index.tar.gz ;;
  mm10)
    wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
    gunzip -c mm10.fa.gz > mm10.full.fa && samtools faidx mm10.full.fa
    samtools faidx mm10.full.fa $(cut -f1 mm10.full.fa.fai | grep -E '^chr([0-9]+|X|Y|M)$') > mm10.fa
    samtools faidx mm10.fa && rm mm10.full.fa*
    bowtie2-build --threads "$THREADS" mm10.fa mm10 && rm mm10.fa.gz ;;
  *) echo "unknown genome $G" >&2; exit 1 ;;
  esac
  echo "$G: fasta $OUT/$G/$G.fa, bowtie2_index $OUT/$G/$G"
done

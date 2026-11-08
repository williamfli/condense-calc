#!/bin/bash
# Download a genome FASTA (+ .fai) and its bowtie2 index into OUTDIR/<genome>/, md5-checked.
#   get_reference.sh OUTDIR THREADS GENOME...      (usually via: condense reference GENOME --refs OUTDIR)
# hg38:      NCBI GRCh38 no-alt analysis set with NCBI's prebuilt bowtie2 index (~7 GB, ~30 min).
# mm10:      ENCODE mm10 no-alt analysis set (66 sequences) with ENCODE's prebuilt bowtie2 index ENCFF309GLL
#            (~4 GB); the same index files as 4DN 4DNFI2493SDN, Park et al.'s "4D_mm10", used if ENCODE fails.
# hg38_ucsc: UCSC hg38.fa as distributed (455 sequences, with alts), Park et al.'s human reference; the index is
#            built here (~1.5 h with 8 threads, ~8 GB of memory).
# Downloads use 16 connections (condense_calc/download.py; one connection to encode-public gives ~0.35 MB/s).
# Runs with macOS's bash 3.2 too.
set -euo pipefail
mkdir -p "$1"; OUT=$(cd "$1" && pwd -P); THREADS=$2; shift 2
PYTHON=${PYTHON:-python}; HERE=$(cd "$(dirname "$0")" && pwd -P)
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
S3=https://encode-public.s3.amazonaws.com
DN=https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/63e22058-79c6-4e24-8231-ca4afac29dda
get() {   # get URL MD5: download into the current folder (resuming a partial file), md5-checked
  "$PYTHON" "$HERE/download.py" "$1" "$2" 16
}
link_index() {   # name the extracted index files <genome>.*.bt2
  for f in $(find . -name '*.bt2' ! -name "$G.*"); do
    ln -sf "$f" $G.$(basename "$f" | grep -oE '(rev\.)?[0-9]\.bt2$'); done
  [ "$(ls $G.*.bt2 | wc -l)" -eq 6 ] || { echo "expected 6 $G index files" >&2; exit 1; }
}
faidx() {   # faidx N: index $G.fa and check that it holds N sequences
  samtools faidx $G.fa
  [ "$(wc -l < $G.fa.fai)" -eq "$1" ] || { echo "$G.fa: expected $1 sequences" >&2; exit 1; }
}
for G in "$@"; do
  mkdir -p "$OUT/$G" && cd "$OUT/$G"
  case $G in
  hg38)
    F=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    get $NCBI/$F.gz a08035b6a6e31780e96a34008ff21bd6
    get $NCBI/$F.bowtie_index.tar.gz 3f493e82d549f1d476491930965c10d0
    gunzip -c $F.gz > hg38.fa && faidx 195
    tar -xzf $F.bowtie_index.tar.gz && link_index
    rm $F.gz $F.bowtie_index.tar.gz ;;
  mm10)
    F=mm10_no_alt_analysis_set_ENCODE.fasta.gz
    get $S3/2016/01/23/b6a24661-be65-48a4-9e4e-44a007f35f5c/$F c3eba224e059be03bdf1de2f30c1a969
    gunzip -c $F > mm10.fa && faidx 66
    { get $S3/2019/10/17/b17e887f-3dad-4380-ab9b-147132feb71e/ENCFF309GLL.tar.gz 76804da5f4ddd74b2eacfe1c2597a6dc &&
      tar -xzf ENCFF309GLL.tar.gz; } ||
    { get $DN/4DNFI2493SDN.bowtie2Index.tar ff76754f3ce07c6a2e8306c7c5119ddd && tar -xf 4DNFI2493SDN.bowtie2Index.tar; }
    link_index
    rm -f $F ENCFF309GLL.tar.gz 4DNFI2493SDN.bowtie2Index.tar ;;
  hg38_ucsc)
    get https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz 1c9dcaddfa41027f17cd8f7a82c7293b
    gunzip -c hg38.fa.gz > $G.fa && faidx 455
    bowtie2-build --threads "$THREADS" $G.fa $G && rm hg38.fa.gz ;;
  *) echo "unknown genome $G" >&2; exit 1 ;;
  esac
  echo "$G: fasta $OUT/$G/$G.fa, bowtie2_index $OUT/$G/$G"
done

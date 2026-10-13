#!/bin/bash
#SBATCH -J condseq_geo_dl
#SBATCH -c 1
#SBATCH -t 0-06:00
#SBATCH -p shared
#SBATCH --mem=1G
#SBATCH -o ../../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../../lab/logs/condseq_pipeline/%x_%j.err
# step 0 (8 files): 1 core, 21 MB peak, 52 min (seff 2026-10-12)
# Download the FASTQ of the steps listed (default: 0 1) from ENA into netscratch, checking md5s.
# Resumable: finished files whose md5 matches are skipped. Submit from condense-seq-pipeline/geo_test/.
#   sbatch download.sbatch 0        # step 0 only (1.28 GB)
#   sbatch download.sbatch 0 1      # steps 0 and 1 (10.1 GB)
# netscratch purges files after 90 days; they can be fetched again with this script. Files get
# the download date (--no-use-server-timestamps), else ENA's old dates trigger the purge at once.
set -uo pipefail
RAW=/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/raw/GSE252941
mkdir -p $RAW
STEPS=" ${@:-0 1} "
FAIL=0
while IFS=$'\t' read -r step run gsm title pairs bytes md5s urls; do
  [ "$step" = step ] && continue
  [[ "$STEPS" == *" $step "* ]] || continue
  IFS=';' read -ra U <<< "$urls"; IFS=';' read -ra M <<< "$md5s"
  for i in "${!U[@]}"; do
    f=$RAW/$(basename "${U[$i]}")
    if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "${M[$i]}" ]; then echo "ok (have) $f"; continue; fi
    wget -q -c --tries=5 --no-use-server-timestamps -O "$f" "https://${U[$i]}" || { echo "DOWNLOAD FAIL $f"; FAIL=1; continue; }
    [ "$(md5sum "$f" | cut -d' ' -f1)" = "${M[$i]}" ] && echo "md5 ok $f" || { echo "MD5 FAIL $f"; FAIL=1; }
  done
done < manifest.tsv
echo "=== done $(date), FAIL=$FAIL, $(du -sh $RAW | cut -f1) in $RAW"
exit $FAIL

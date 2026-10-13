#!/bin/bash
#SBATCH -J condseq_geo_dl
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p shared
#SBATCH --mem=1G
#SBATCH -o ../../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../../lab/logs/condseq_pipeline/%x_%j.err
# step 0 (8 files, one at a time): 1 core, 21 MB peak, 52 min (seff 2026-10-12). ENA gave each
# connection only ~0.4 MB/s, so files are now fetched 4 at a time (wget waits on the network, so
# one core is enough).
# Download the FASTQ of the steps listed (default: 0 1) from ENA into netscratch, checking md5s.
# Files whose md5 already matches are skipped. ENA sometimes serves a file at full length with a few
# hundred bytes wrong (2026-10-13: 14 of 40 step-1 files; one file was bad in all 4 downloads,
# at different offsets each time, also when written to node-local /scratch), so each file is
# downloaded again, up to 5 times, until its md5 matches. Submit from condense-seq-pipeline/geo_test/.
#   sbatch download.sbatch 0        # step 0 only (1.28 GB)
#   sbatch download.sbatch 0 1      # steps 0 and 1 (10.1 GB)
# netscratch purges files after 90 days; they can be fetched again with this script. Files get
# the download date (--no-use-server-timestamps), else ENA's old dates trigger the purge at once.
set -uo pipefail
RAW=/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/raw/GSE252941
PARALLEL=4
mkdir -p $RAW
STEPS=" ${@:-0 1} "

fetch() {   # fetch FILE MD5 URL; up to 5 attempts, a fresh download each time
  local f=$1 md5=$2 url=$3 try
  if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then echo "ok (have) $f"; return 0; fi
  for try in 1 2 3 4 5; do
    rm -f "$f"
    wget -q --tries=5 --no-use-server-timestamps -O "$f" "https://$url" || { echo "download error (try $try) $f"; continue; }
    if [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then echo "md5 ok (try $try) $f"; return 0; fi
    echo "md5 mismatch (try $try) $f"
  done
  rm -f "$f"; echo "FAIL after 5 tries $f"; return 1
}
export -f fetch

while IFS=$'\t' read -r step run gsm title pairs bytes md5s urls; do
  [ "$step" = step ] && continue
  [[ "$STEPS" == *" $step "* ]] || continue
  IFS=';' read -ra U <<< "$urls"; IFS=';' read -ra M <<< "$md5s"
  for i in "${!U[@]}"; do
    printf '%s\t%s\t%s\n' "$RAW/$(basename "${U[$i]}")" "${M[$i]}" "${U[$i]}"
  done
done < manifest.tsv | xargs -P $PARALLEL -L 1 -d '\n' bash -c 'IFS=$'"'"'\t'"'"' read -r f m u <<< "$0"; fetch "$f" "$m" "$u"'
FAIL=$?
echo "=== done $(date), exit=$FAIL (non-zero: some file failed), $(du -sh $RAW | cut -f1) in $RAW"
exit $FAIL

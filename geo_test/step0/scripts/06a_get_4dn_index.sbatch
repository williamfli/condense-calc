#!/bin/bash
#SBATCH -J condseq_s0_4dnidx
#SBATCH -c 1
#SBATCH -t 0-03:00
#SBATCH -p shared
#SBATCH --mem=1G
#SBATCH -o ../../../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../../../lab/logs/condseq_pipeline/%x_%j.err
# Step 0, mechanism test 3: fetch the 4DN GRCh38 bowtie2 index (4DNFIMQPTYDY, 4.17 GB; listed as released),
# the reference Park's GEO run names as "4D_hg38". Into netscratch (throwaway), md5-checked against
# the 4DN portal's record. Submit from condense-seq-pipeline/geo_test/step0/.
# NOT RUN: 4DN answers 403 to anonymous downloads (checked 2026-10-13), even for released reference
# files. Needs a 4DN access key (portal profile -> Access Keys), passed through the environment:
#   FOURDN_KEY=... FOURDN_SECRET=... sbatch --export=ALL scripts/06a_get_4dn_index.sbatch
set -euo pipefail
: "${FOURDN_KEY:?set FOURDN_KEY and FOURDN_SECRET (4DN access key)}" "${FOURDN_SECRET:?}"

D=/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/ref_4dn/hg38
mkdir -p $D && cd $D
wget -q -c --no-use-server-timestamps --user "$FOURDN_KEY" --password "$FOURDN_SECRET" \
  https://data.4dnucleome.org/files-reference/4DNFIMQPTYDY/@@download/4DNFIMQPTYDY.bowtie2Index.tar
[ "$(md5sum 4DNFIMQPTYDY.bowtie2Index.tar | cut -d' ' -f1)" = 6a30578bdc3218b1caa1ac44809c4117 ] || { echo MD5 FAIL; exit 1; }
tar -xmf 4DNFIMQPTYDY.bowtie2Index.tar && rm 4DNFIMQPTYDY.bowtie2Index.tar
find $D -name '*.bt2' -exec ls -la {} \;

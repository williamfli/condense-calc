#!/bin/bash
#SBATCH -J condseq_equiv
#SBATCH -c 2
#SBATCH -t 0-04:00
#SBATCH -p test
#SBATCH --mem=24G
#SBATCH -o ../lab/logs/condseq_pipeline/%x_%j.out
#SBATCH -e ../lab/logs/condseq_pipeline/%x_%j.err
# Run Park's original Python 2 scripts (external/condense-seq/prepro_scripts) on the synthetic
# BAMs made by tests/run_synthetic.sbatch, with the GEO run's arguments, then compare with ours.
# Submit from condense-seq-pipeline/ after the synthetic run has finished.
set -euo pipefail
module load python/3.12.11-fasrc02
mamba activate condseq                       # samtools, and python3 for the comparison
module load python/2.7.16-fasrc01            # python2 for the original scripts
ORIG=$(realpath ../external/condense-seq/prepro_scripts)
PIPE=$(pwd)
OUT=$(realpath ../lab/outputs)/condense-seq-pipeline/synthetic
REF=$OUT/input/ref/sim.fa
CHR="chr1 chr2 chrX"
W=$OUT/equiv/original; mkdir -p $W; cd $W
PY2=$(command -v python2)
PY3=/n/holylabs/zhuang_lab/Lab/wfli/conda/envs/condseq/bin/python

# Park-format titration file (10 tab-separated columns: conc .. survival fraction (8th) .. tnum)
awk 'NR==1{print "conc\tc1\tc2\tc3\tc4\tc5\tc6\tfrac\tc8\ttnum"; next}
     {printf "%s\t0\t0\t0\t0\t0\t0\t%s\t0\t%s\n", $2, $3, $1}' $OUT/input/titration/SIM_NCP_sp.tsv > tit.csv

BAMS=$(ls $OUT/run/bam/*.bam)
for b in $BAMS; do n=$(basename $b .bam)
  for w in 1000 10000; do kb=$((w/1000))kb
    $PY2 $ORIG/bincount_ver3.py $b $REF -o ${n}_${kb} -w $w --chr $CHR --gc --tlen 2>/dev/null
  done
done
for b in $(ls $OUT/run/bam/*_deep.bam); do n=$(basename $b .bam)
  for c in $CHR; do
    $PY2 $ORIG/coverage_ver3.py $b $REF --chr $c --skip -o ${n}_${c} 2>/dev/null
  done
done
INP=SIM_NCP_sp_0_1rep_deep
for c in $CHR; do
  $PY2 $ORIG/NCPpeak_ver3.py $OUT/run/bam/$INP.bam --chr $c --skip -o ${INP}_${c} 2>/dev/null
done
for b in $(ls $OUT/run/bam/*_deep.bam); do n=$(basename $b .bam)
  for c in $CHR; do
    $PY2 $ORIG/NCPcov_ver3.py ${INP}_${c}_peak.gtab.gz ${n}_${c}_cov.gtab.gz --Nlen 171 --peak-choice input --chr $c -o ${n}_${c} 2>/dev/null
    $PY2 $ORIG/bindata_ver3.py ${n}_${c}_cov.gtab.gz $REF --chr $c --Bsize 171 --Bstep 25 -o ${n}_${c}_171win25step 2>/dev/null
  done
done
# scores: binned (shallow and deep) and per nucleosome
for b in $BAMS; do n=$(basename $b .bam); t=$(echo $n | cut -d_ -f4); [ "$t" = 0 ] && continue
  inp=$(echo $n | sed "s/_sp_${t}_/_sp_0_/")
  for kb in 1kb 10kb; do
    $PY2 $ORIG/NCPscore_ver4.py ${n}_${kb}_bin.gtab.gz -i ${inp}_${kb}_bin.gtab.gz -t tit.csv --tnum $t --chr $CHR --numc 2>/dev/null
  done
  case $n in *_deep)
    $PY2 $ORIG/NCPscore_ver4.py $(for c in $CHR; do echo ${n}_${c}_Ncov.gtab.gz; done) \
      -i $(for c in $CHR; do echo ${inp}_${c}_Ncov.gtab.gz; done) -t tit.csv --tnum $t --chr $CHR --numc 2>/dev/null ;;
  esac
done
# nucleosome numbers (10 kb, shallow series) -> concatenate -> c-half fit (logistic_fit_ver3.py is Python 3)
for t in 0 1 2 3 4 5; do
  $PY2 $ORIG/NCPnum_ver2.py SIM_NCP_sp_${t}_1rep_10kb_bin.gtab.gz -t tit.csv --tnum $t --mscale 1 --chr $CHR 2>/dev/null
done
$PY2 $ORIG/concat_gtab.py $(for t in 0 1 2 3 4 5; do echo SIM_NCP_sp_${t}_1rep_10kb_num.gtab.gz; done) \
  -c 3 3 3 3 3 3 -o SIM_NCP_sp_1rep_10kb_num 2>/dev/null
$PY3 $ORIG/logistic_fit_ver3.py SIM_NCP_sp_1rep_10kb_num.gtab.gz -t tit.csv --tnum 0 1 2 3 4 5 --model sigmoid \
  --min_rsq 0.5 --min_top 1 --max_top 1 --min_bottom 0 --max_bottom 0 --min_rate 0 --max_rate 100 -o SIM_NCP_sp_1rep_10kb
ls | wc -l

# our coverage files, for the coverage comparison (the workflow keeps coverage in memory)
mkdir -p $OUT/equiv/ours && cd $PIPE
for b in $(ls $OUT/run/bam/*_deep.bam); do n=$(basename $b .bam)
  for c in $CHR; do
    PYTHONPATH=$PIPE $PY3 -m condseq coverage --bam $b --fasta $REF --chrom $c -o $OUT/equiv/ours/${n}_${c}_cov.gtab.gz
  done
done
PYTHONPATH=$PIPE $PY3 tests/equivalence/compare.py $OUT

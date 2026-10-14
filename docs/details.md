# Processing details

condseq is a Python 3 port of the preprocessing scripts of Park et al. 2025
([spark159/condense-seq](https://github.com/spark159/condense-seq) `prepro_scripts`, commit `a8d2dd8`),
with the settings of the script that made the GEO tables (`submit_scripts/batch_submit_ver7.py`).

## Steps

`<name>` is `<cell>_<sample>_<agent>_<tnum>_<rep>rep[_deep]`, `<series>` the same without `<tnum>`.

| step | original script | settings | output |
|---|---|---|---|
| align | `bowtie-submit` | bowtie2 end-to-end, paired, `--trim-to 50` | `bam/<name>.bam` |
| bincount | `bincount_ver3.py` | 1 kb and 10 kb bins, fragment centres, no length filter | `bin/` |
| score_bin | `NCPscore_ver4.py --numc` | each point vs its input | `score/` |
| table_bin | `concat_gtab.py` | one column per point | `tables/<series>_<kb>_score.gtab.gz` |
| peaks | `NCPpeak_ver3.py` | deep input only; 126–184 bp fragments; 147 bp nucleosomes overlapping ≤ 40 bp | `peaks/` |
| areas | `coverage_ver3.py`, `NCPcov_ver3.py`, `bindata_ver3.py` | coverage in 171 bp at each dyad; 171 bp windows every 25 bp | `areas/` |
| score_ncov, table_ncov | `NCPscore_ver4.py`, `concat_gtab.py` | normalized over all chromosomes | `tables/<series>_score/` |
| num_bin, table_num | `NCPnum_ver2.py` | 1.6 × 10¹² input nucleosomes × agent scale | `tables/<series>_10kb_num.gtab.gz` |
| fit | `logistic_fit_ver3.py` | 10 kb, series with ≥ 5 points; sigmoid, top 1, bottom 0, rate 0–100, differential evolution (seed 3); R² ≥ 0.5 | `tables/<series>_10kb_4PL.gtab.gz` |

Score: −ln[(s+1)/Σ(s+1) ÷ (c+1)/Σ(c+1)] − ln *f*, summed over rows with input *c* > 0; rows without
input reads are dropped.

Reproduced as in the original:
- **Reads:** each fragment is counted once, from its forward read. Kept: proper pairs with the mate
  mapped and edit distance ≤ 10. No MAPQ filter and no duplicate removal.
- **Fragment centre:** even lengths go to the left centre base in `bincount`; the peak caller splits
  them 0.5/0.5.
- **Peak caller:** greedy and unsmoothed. Equal dyad counts are ordered as Python 2 sorted them.
  `condseq peaks --fast` uses a plain tie rule instead.
- **Numbers in text output:** rounded and printed as Python 2 did, so files compare byte for byte.
- **Not available:** `bowtie-submit_ver2`, which the GEO driver calls, is missing from the authors'
  repository; the align step uses `bowtie-submit` with the driver's `trim_to = 50`.

Sliding-window areas (`areas/*_Bdata*`) are computed but not scored; GEO has no table of them.
`condseq score` can score them.

## Validation

**Against the original scripts** (`tests/equivalence/run.sh`). A simulated experiment of 600 kb is
used: 3 chromosomes, a duplicated 2 kb segment, a 5 kb N gap, 3,127 nucleosomes, a 6-point shallow
series and a 3-point deep series. Every table is identical:

| step | files | identical |
|---|---|---|
| bincount | 18 | 17,820/17,820 cells |
| score (bins) | 14 | 4,585/4,585 |
| coverage | 9 | 1,758,580/1,758,580 |
| ncov | 9 | 9,381/9,381 |
| bindata | 9 | 71,307/71,307 |
| score (nucleosomes) | 6 | 6,254/6,254 |
| peaks | 3 | 3,127/3,127 |
| num | 6 | 360/360 |
| fit | 1 | 420/420 |

**Against the simulated truth** (`tests/synthetic/check_truth.py`): 10 kb scores r = 0.98–0.998
(points 2–5); 1 kb deep r = 0.96–0.99; per-nucleosome r = 0.89–0.95, with 3,118 of 3,127 dyads
found within 20 bp; c½ r = 0.991.

**Against GEO (real reads).** We re-ran GEO samples from their raw reads (hg38, NCBI no-alt
analysis set):
- H1 PEG, replicate 2: r = 0.964 with the GEO 10 kb table; r = 0.995 in bins without multi-mapping
  reads.
- H1 spermine, replicate 2, nine points: r = 0.94–0.98, with no offset at any point.

The remaining differences are in repeats. A multi-mapping read's position depends on the genome
index, and Park et al. used the 4DN hg38 index. `reads: {min_mapq: 10}` removes these reads, at the
cost of no longer matching GEO.

## Resources

Measured on hg38:

| step | memory | time |
|---|---|---|
| align | ~6 GB | about 0.8 M read pairs per minute with 8 threads |
| bincount | ~3 GB plus 0.6 × the BAM size | |
| score_bin | 4.4 GB at 1 kb | |
| fit | 0.5 GB | ~2.5 h for 10 points (single core) |

All other steps need ≤ 2 GB. Memory requests double on the retry after an out-of-memory failure.
Deep samples (~300 M pairs) take most of a day to align.

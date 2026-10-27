# Processing details

condense-calc is a Python 3 port of the preprocessing scripts of Park et al. 2025
([spark159/condense-seq](https://github.com/spark159/condense-seq) `prepro_scripts`, commit `a8d2dd8`),
with the settings of the script that made the GEO tables (`submit_scripts/batch_submit_ver7.py`).
The defaults change a few of them (`condense_calc/settings.py`); `--park` (`park: true`) restores all.

## Steps

`<name>` is one library: `<series>_<point>`, where the point numbers the concentrations of a series
from 0 (the input). A sheet can set it in a `point` column. With a config (below) the names follow
GEO: `<cell>_<sample>_<agent>_<tnum>_<rep>rep[_deep]`.

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
input reads are written as NA (`--park`: dropped).

Reproduced as in the original, with the default's differences:
- **Reads:** each fragment is counted once, from its forward read. Kept: proper pairs with the mate
  mapped and edit distance ≤ 10, MAPQ ≥ 10 (`--park`: no MAPQ filter). No duplicate removal.
- **Fragment centre:** even lengths go to the left centre base in `bincount`; the peak caller splits
  them 0.5/0.5. The default takes the centre from the fragment's interval; `--park` keeps the
  original's formula, one base off for odd fragments counted from the right-hand read (never the case
  with bowtie2 end-to-end pairs), and its wrap of a centre before position 0 to the last bin.
- **Peak caller:** greedy and unsmoothed. Equal dyad counts are taken from the lowest position up;
  `--park` orders them as Python 2 sorted them.
- **Numbers in text output:** rounded and printed as Python 2 did, so files compare byte for byte.
- **Not available:** `bowtie-submit_ver2`, which the GEO driver calls, is missing from the authors'
  repository; the align step uses `bowtie-submit` with the driver's `trim_to = 50`.

Sliding-window areas (`areas/*_Bdata*`) are computed only with `bdata: true` and not scored; GEO has
no table of them. `condense-calc score` can score them.

## Validation

**Against the original scripts** (`tests/equivalence/run.sh`, in Park mode). A simulated experiment of 600 kb is
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
(points 2–5); 1 kb deep r = 0.96–0.99; per-nucleosome r = 0.89–0.95, with 3,109 dyads called, all
within 20 bp of a true one (the 18 others lie in the duplicated segment, where MAPQ ≥ 10 leaves no
reads; with `--park` 3,118 of 3,127 called dyads are within 20 bp); c½ r = 0.991.

**Against GEO (real reads).** We re-ran GEO samples from their raw reads (hg38, NCBI no-alt
analysis set):
- H1 PEG, replicate 2: r = 0.964 with the GEO 10 kb table; r = 0.995 in bins without multi-mapping
  reads.
- H1 spermine, replicate 2, nine points: r = 0.94–0.98, with no offset at any point.

These runs used Park's settings with ENA reads and the NCBI index. The remaining differences are in
repeats, where reads map to several places: where bowtie2 places them depends on the reference, the
read names and the base qualities. With Park's inputs (reads as `condense fetch` makes them, Park's
references and bowtie2 2.4.x, as `--park` sets, and full-precision titrations) the H1 and E14 tables
tested were reproduced byte for byte. The default MAPQ filter removes such reads.

## Resources

Measured on hg38:

| step | memory | time |
|---|---|---|
| align | ~6.3 GB | ~1 M read pairs per minute with 8 threads; a deep library (~330 M pairs) ~2 h with 24 |
| bincount | ~3 GB plus 0.6 × the BAM size (deep: 8.4 GB used) | deep: ~20 min |
| score_bin | 4.4–4.9 GB at 1 kb | ~1 min |
| peaks (per chromosome) | ≤ 11.5 GB (deep, chr1) | ≤ 5 min |
| areas (per chromosome) | ≤ 6.5 GB (deep) | ≤ 3 min |
| score_ncov | ≤ 14.8 GB (deep) | ≤ 7 min |
| fit | 0.5 GB | ~2.5 h for 10 points (single core) |

Deep figures: H1, 3 libraries of ~320 M pairs each, 2026-10-19. Other steps need ≤ 2 GB. Memory requests
double on the retry after an out-of-memory failure, and time limits double on the retry after a timeout
(up to `max_runtime`, 3 days). A deep series with per-nucleosome scores takes ~3–4 h from FASTQ.

## All settings

Given with `--settings FILE` (YAML), or in a config file.

| key | default | |
|---|---|---|
| `threads` | 8 | bowtie2 threads per library |
| `bin_sizes` | `[1000, 10000]` | |
| `park` | `false` | `true`: Park's processing (`--park`); the keys below then take Park's values unless given |
| `reads: {min_mapq: 10}` | MAPQ ≥ 10 (Park: 0, no filter) | |
| `reads: {dedup: false}` | duplicates kept, as Park | |
| `reads: {mm: 10}` | | maximum edit distance |
| `score: {missing_as_nan: true}` | bins without input reads are NA rows (Park: dropped) | |
| `fastq_check` | `lite` (Park: `park`) | `off`, `lite` (stop on SRA Lite), `park` (also read names) |
| `bowtie2` | `bowtie2`: this environment's (Park: `2.4.1`) | a version (installed once under the refs folder) or a path |
| `bdata` | `false` | `true`: also sliding-window areas of deep series |
| `mscale` | `{Ca: 5}`, others 1 | molecule-number scale per agent |
| `fit` | `{bin_size: 10000, min_tnum: 5, model: sigmoid, method: evolution, min_rsq: 0.5}` | c½ fit |
| `trim_to` | 50 | reads are trimmed to 50 bp before alignment |
| `qc` | `default` | `off`, `default`, `verbose` |
| `qc_thresholds` | | overrides, e.g. `{min_replicate_rho: 0.7}`; the thresholds are listed in the report |
| `chroms` | hg38/mm10-like: 1–N, X, Y; other genomes: every sequence ≥ 1% of the longest, except the mitochondrial (Park: no chrY for GM and mCD8T) | chromosomes to score; e.g. without chrY for female cells |
| `max_runtime` | 4320 | minutes; the longest time limit a retried step may ask for |

## Config files (several genomes, GEO names)

`condense config.yaml` runs from a config instead of a sample sheet. Use this for libraries from
several genomes, or to get file names as in GEO. Relative paths are read from the config's folder.

```yaml
samples: samples.tsv        # cell, sample, agent, tnum, rep, depth, fastq1, fastq2 (lanes comma-separated)
titration_dir: titration    # <cell>_<sample>_<agent>.tsv: tnum, conc, survival_fraction
genomes:                    # `condense reference hg38` prints these paths
  hg38: {fasta: refs/hg38/hg38.fa, bowtie2_index: refs/hg38/hg38}
cells:
  H1: hg38                  # or {genome: hg38, chroms: [chr1, chr2]}
```

`depth: deep` gives per-nucleosome scores for that series. Any setting above can be added.

## Tests

```bash
pytest tests/test_units.py
bash tests/synthetic/run.sh /tmp/cc-test             # simulated experiment, ~2 min
python tests/synthetic/check_truth.py /tmp/cc-test   # recovery of the simulated scores and c-half
bash tests/synthetic/run_sheet.sh /tmp/cc-test       # same run from a sample sheet: identical tables
bash tests/equivalence/run.sh /tmp/cc-park PARK_REPO PYTHON2   # --park run: identical to the original scripts
```

# Settings, config files and steps

condense-calc is a Python 3 port of the preprocessing scripts of Park & Ha, Nature, 2025
([spark159/condense-seq](https://github.com/spark159/condense-seq) `prepro_scripts`, commit `a8d2dd8`),
with the settings of the script that made the GEO tables (`submit_scripts/batch_submit_ver7.py`).
Where the defaults differ from Park et al.'s processing, and what `--park` restores: [park.md](park.md).

## All settings

In a YAML file given with `--settings FILE` (sample-sheet run), or in the config file (config run).
`--park`, `--qc` and (sample-sheet runs) `--nucleosomes` set `park`, `qc` and `nucleosomes`. Settings you give
override `--park`. `--settings`, `--nucleosomes` and `-g` have no effect on a config run.

| setting | default | `--park` | |
|---|---|---|---|
| `park` | `false` | `true` | Park et al.'s processing ([park.md](park.md)) |
| `reads: {min_mapq}` | 10 | 0 | minimum MAPQ of a fragment's forward-strand read; 0: no filter |
| `reads: {mm}` | 10 | | maximum edit distance (NM) |
| `reads: {dedup}` | `false` | | `true`: drop fragments with the same chromosome, start and length |
| `score: {missing_as_nan}` | `true` | `false` | rows without input reads: NA (`true`) or dropped (`false`) |
| `fastq_check` | `lite` | `park` | `off`; `lite`: stop on SRA Lite reads; `park`: also on read names not in the instrument's form |
| `bowtie2` | `bowtie2` | `2.4.1` | `bowtie2`: the environment's (2.5.5); a full version (x.y.z): installed once in `<refs>/envs/`; or the path of a bowtie2 |
| `bdata` | `false` | | `true`: deep series also get sliding-window areas (not scored) |
| `threads` | 8 | | bowtie2 threads per library |
| `trim_to` | 50 | | reads are trimmed to this length before alignment |
| `bin_sizes` | `[1000, 10000]` | | bin sizes in bp |
| `mscale` | `{Ca: 5}` | | nucleosome-number scale per agent; others 1 |
| `fit` | `{bin_size: 10000, min_tnum: 5, model: sigmoid, method: evolution, min_rsq: 0.5}` | | c½ fit; `min_tnum` counts the input; `model`: `sigmoid` or `hill`; `method`: `evolution` or `curve_fit` |
| `qc` | `default` | | `off`, `default`, `verbose` (adds duplicates, GC, a score profile, dinucleotide periodicity) |
| `qc_thresholds` | `{}` | | overrides, e.g. `{min_replicate_rho: 0.7}`; the report lists every threshold |
| `max_runtime` | 4320 | | minutes; the longest time limit a retried step may ask for (Slurm) |
| `chroms` | see below | GM, mCD8T: no chrY | chromosomes to score; in a config file, per cell under `cells:` |
| `nucleosomes` | `false` | | sample-sheet run: per-nucleosome scores for every series (config run: `depth: deep`) |
| `refs` | `$CONDENSE_REFS`, else `~/.cache/condense-calc` | | where genomes and bowtie2 versions are kept (`--refs`) |

Default `chroms`: for human, mouse and genomes like them (chromosomes named 1..N, X, Y), exactly those, in
that order; for other genomes, every sequence at least 1% as long as the longest, except the mitochondrial one.
Config files also take `samples`, `titration_dir`, `genomes` and `cells` (below). `sheet` and `genome` are
written by `condense` for a sample-sheet run; `align_threads` is an older name of `threads`.

Example, used as `condense samples.tsv --settings settings.yaml`:

```yaml
reads: {min_mapq: 0}                # keep multi-mapping reads
chroms: [chr1, chr2, chr3, chrX]    # score only these, e.g. to leave out chrY
```

## Config files

`condense config.yaml` runs from a config file instead of a sample sheet. Use it for libraries of several
genomes, or for file and column names as on GEO. Relative paths are read from the config's folder.
Example: H1 PEG replicate 2 of GSE252941.

```yaml
samples: libraries.tsv
titration_dir: titration
genomes:                    # `condense reference hg38_ucsc` prints these paths
  hg38: {fasta: ~/.cache/condense-calc/hg38_ucsc/hg38_ucsc.fa, bowtie2_index: ~/.cache/condense-calc/hg38_ucsc/hg38_ucsc}
cells:
  H1: hg38                  # or {genome: hg38, chroms: [chr1, chr2]}
```

`libraries.tsv`, tab-separated (spaces do not work), one row per library (lanes comma-separated;
`depth: deep` adds per-nucleosome scores):

```
cell  sample  agent  tnum  rep  depth    fastq1                                                     fastq2
H1    NCP     PEG    0     2    shallow  fastq/SRR31626988_1.fastq.gz,fastq/SRR31626989_1.fastq.gz  fastq/SRR31626988_2.fastq.gz,fastq/SRR31626989_2.fastq.gz
H1    NCP     PEG    6     2    shallow  fastq/SRR31626986_1.fastq.gz,fastq/SRR31626987_1.fastq.gz  fastq/SRR31626986_2.fastq.gz,fastq/SRR31626987_2.fastq.gz
```

`titration/<cell>_<sample>_<agent>.tsv`, here `titration/H1_NCP_PEG.tsv`:

```
tnum  conc  survival_fraction
0     0     1
6     9.5   0.097744361
```

`tnum` numbers the concentrations as on GEO (0: the input); `survival_fraction` is the soluble fraction.
This gives `results/tables/H1_NCP_PEG_2rep_10kb_score.gtab.gz` with the column `H1_NCP_PEG_6_2rep`, as on GEO.
Any setting above can be added to the config.

## Steps

`<name>` is one library: `<series>_<point>`, where the point numbers the concentrations of a series from 0
(the input); a sheet can set it in a `point` column. With a `rep` column the series is `<series>_rep<N>`.
With a config file the names follow GEO: `<cell>_<sample>_<agent>_<point>_<rep>rep[_deep]`.

| step | Park's script | what it does | output |
|---|---|---|---|
| align | `bowtie-submit` | FASTQ check (`fastq_check`), then bowtie2 end-to-end, paired, `--trim-to 50` | `bam/<name>.bam` |
| bincount | `bincount_ver3.py` | fragments per bin (`bin_sizes`), at the fragment centre, no length filter | `bin/` |
| score_bin | `NCPscore_ver4.py --numc` | each point against its input | `score/` |
| table_bin | `concat_gtab.py` | one column per point | `tables/<series>_<kb>_score.gtab.gz` |
| peaks | `NCPpeak_ver3.py` | deep series, on the input: dyads of 126–184 bp fragments; 147 bp nucleosomes overlapping ≤ 40 bp | `peaks/` |
| areas | `coverage_ver3.py`, `NCPcov_ver3.py`; `bindata_ver3.py` | coverage in 171 bp at each dyad; with `bdata: true` also in 171 bp windows every 25 bp | `areas/` |
| score_ncov, table_ncov | `NCPscore_ver4.py`, `concat_gtab.py` | per nucleosome, normalized over all chromosomes | `tables/<series>_score/` |
| num_bin, table_num | `NCPnum_ver2.py` | fitted series: nucleosome numbers, 1.6 × 10¹² input nucleosomes × `mscale` | `tables/<series>_10kb_num.gtab.gz` |
| fit | `logistic_fit_ver3.py` | c½ at 10 kb, series with ≥ 5 points; sigmoid, top 1, bottom 0, rate 0–100, differential evolution (seed 3); R² ≥ 0.5 | `tables/<series>_10kb_4PL.gtab.gz` |
| QC, run record | | QC report and flags; settings, references and versions | `qc/`, `run_info.json` |

- **Score:** −ln[(s+1)/Σ(s+1) ÷ (c+1)/Σ(c+1)] − ln *f*, for supernatant count *s*, input count *c* and soluble
  fraction *f*; the sums run over rows with *c* > 0. Rows without input reads are NA (`--park`: dropped).
- **Reads:** each fragment is counted once, from its forward-strand read. Kept: proper pairs with the mate
  mapped, edit distance (NM) ≤ 10 and MAPQ ≥ `min_mapq`. No duplicate removal.
- **Fragment centre:** even lengths go to the left of the two centre bases in `bincount`; the peak caller
  splits them 0.5/0.5. The peak caller is greedy and does no smoothing.
- **Text output:** numbers are rounded and printed as Python 2 did, so files compare byte for byte.
- `bowtie-submit_ver2`, which Park's driver calls, is missing from Park's repository; the align step uses
  `bowtie-submit` with the driver's `trim_to = 50`.

Not carried over from Park's scripts: per-base coverage files and the full dyad-signal file (both computed in
memory), scores of sliding windows, the log2/log10 score options, and Park's downstream scripts
(fragment-length plots, occupancy, motifs, profiles, annotation).

## Resources

Measured on hg38. Deep: libraries of ~320 M read pairs.

| step | memory | time |
|---|---|---|
| align | ~6.3 GB | ~0.8–1 M read pairs per minute with 8 threads; a deep library ~2 h with 24 |
| bincount | ~3 GB plus 0.6 × the BAM size (deep: 8.4 GB) | deep: ~20 min |
| score_bin | 4.4–4.9 GB at 1 kb | ~1 min |
| peaks (per chromosome) | ≤ 11.5 GB (deep, chr1) | ≤ 5 min |
| areas (per chromosome) | ≤ 6.5 GB (deep) | ≤ 3 min |
| score_ncov | ≤ 14.8 GB (deep) | ≤ 7 min |
| fit | 0.5 GB | ~2.5 h for 10 points (one core) |

Other steps need ≤ 2 GB. With `--slurm`, a failed step is retried once with twice the memory and twice the time
limit (time capped at `max_runtime`). A deep series with per-nucleosome scores takes ~3–4 h from FASTQ.
Genomes, once: hg38 ~7 GB download (~30 min); `hg38_ucsc` builds its index (~1.5 h with 8 threads, ~8 GB).

## Tests

From the repo root, in the condense-calc environment:

```bash
pytest tests/test_units.py
bash tests/synthetic/run.sh /tmp/cc-test             # simulated experiment, ~1–2 min
python tests/synthetic/check_truth.py /tmp/cc-test   # recovery of the simulated scores and c½
bash tests/synthetic/run_sheet.sh /tmp/cc-test       # the same run from a sample sheet: the same values
bash tests/equivalence/run.sh /tmp/cc-park PARK_REPO PYTHON2   # Park mode against Park's scripts
```

`PARK_REPO` is a checkout of spark159/condense-seq at `a8d2dd8`; `PYTHON2` defaults to `python2`.
The equivalence test runs condense-calc in Park mode on a simulated 600 kb experiment (a shallow and a deep
series), then Park's Python 2 scripts on the same BAMs: every intermediate and final table is identical.
`check_truth.py` compares the tables with the simulated truth (10 kb scores r = 0.98–0.998, per nucleosome
r = 0.89–0.95, c½ r = 0.99). Against GEO: [park.md](park.md#reproducing-the-geo-tables-exactly).

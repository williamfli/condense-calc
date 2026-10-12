# condense-seq-pipeline

From condense-seq raw reads (paired FASTQ) to condensability tables in the GEO GSE252941 format:
binned scores (1 kb, 10 kb), per-nucleosome scores, nucleosome numbers and c½ fits.

It is a Python 3 port of the Park 2025 preprocessing (`../external/condense-seq/prepro_scripts`,
commit `a8d2dd8`). The settings come from the script that made the GEO files
(`submit_scripts/batch_submit_ver7.py` and the `*-submit` job scripts). **On simulated reads,
every stage gives byte-identical output to Park's original scripts** (see Validation).
Departures from Park are opt-in options, off by default.

The polymer simulation (Park Supp. Note 2) is out of scope here.

Status (2026-10-07): built and validated on simulated reads only. **Not yet run on real reads;
no reference genomes are installed yet** (see "Running on real data").

## Layout

```
condense-seq-pipeline/
├── condseq/                  Python package, one module per group of stages
│   ├── fragments.py          BAM -> fragments with Park's read filters
│   ├── counts.py             bincount, coverage, peaks, ncpcov, bindata
│   ├── score.py              score, numbers, concatenate
│   ├── fit.py                c½ logistic fit
│   ├── titration.py          titration table reader (ours, or Park's 10-column format)
│   ├── gtab.py               the gtab table format, chromosome order, Python 2 number formatting
│   └── cli.py                `python -m condseq <stage>`
├── workflow/Snakefile        sample sheet -> GEO-format tables
├── workflow/profiles/slurm/  Snakemake profile: one Slurm job per step
├── config/config.example.yaml
├── env/                      conda env `condseq` (spec + lock)
├── slurm/                    env build, reference setup
├── tools/                    titration tables from Park's Source Data
├── geo_test/                 first real-data check against GEO (prepared, not run)
└── tests/
    ├── test_units.py         small hand-checked cases (pytest)
    ├── simulate.py           synthetic experiment with a known answer
    ├── run_synthetic.sbatch  simulate + run the whole workflow
    ├── check_truth.py        compare outputs with the simulated truth
    └── equivalence/          run Park's Python 2 scripts on the same BAMs and compare
```
Everything large goes to lab storage:
- outputs: `../lab/outputs/condense-seq-pipeline/<run>/`
- references: `../lab/data/ref/`
- logs: `../lab/logs/condseq_pipeline/`

## Stages

`<name>` is `<cell>_<sample>_<agent>_<tnum>_<rep>rep[_deep]`, the GEO naming; `tnum` 0 is the input.

| Step (rule) | Port of | GEO-run settings | Output |
|---|---|---|---|
| align | `bowtie-submit` | bowtie2 end-to-end, `--trim-to 50`, paired; then sorted + indexed | `bam/<name>.bam` |
| bincount | `bincount_ver3.py` | `-w 1000` / `-w 10000 --gc --tlen`, no length filter | `bin/<name>_<kb>_bin.gtab.gz` |
| score_bin | `NCPscore_ver4.py --numc` | each titration point vs its input | `score/<name>_<kb>_score.gtab.gz` |
| table_bin | `concat_gtab.py` | one column per titration point | **`tables/<series>_<kb>_score.gtab.gz`** |
| peaks | `NCPpeak_ver3.py` | input only, deep only, 126–184 bp fragments, 147 bp nucleosome, ≤ 40 bp overlap | `peaks/<input>_<chr>_peak.gtab.gz` |
| areas | `coverage_ver3.py` + `NCPcov_ver3.py` + `bindata_ver3.py` | all fragment lengths; 171 bp per dyad (`--Nlen 171`); windows 171 bp / 25 bp | `areas/<name>_<chr>_Ncov…`, `…_171win25step_Bdata…` |
| score_ncov | `NCPscore_ver4.py --numc` | normalized over all chromosomes | `score/<name>_<chr>_score.gtab.gz` |
| table_ncov | `concat_gtab.py` | per chromosome (GEO's `*_deep_score.tar.gz` content) | **`tables/<series>_score/<series>_<chr>_score.gtab.gz`** |
| num_bin, table_num | `NCPnum_ver2.py`, `concat_gtab.py` | 1.6×10¹² input nucleosomes × agent scale (Ca: 5) | `tables/<series>_10kb_num.gtab.gz` |
| fit | `logistic_fit_ver3.py` | 10 kb, series with ≥ 5 points; sigmoid; differential evolution (seed 3); top = 1, bottom = 0, rate 0–100; R² ≥ 0.5 | **`tables/<series>_10kb_4PL.gtab.gz`** |

`<series>` is `<cell>_<sample>_<agent>_<rep>rep[_deep]`.

**The score**, per bin or dyad, at a titration point with soluble fraction *f*:
−ln[(s+1)/Σ(s+1) ÷ (c+1)/Σ(c+1)] − ln f.
- *s* is the supernatant count and *c* the input count.
- The sums run over rows with *c* > 0.
- Rows whose input is 0 are dropped.

**Details reproduced on purpose:**
- **Reads:**
  - Each fragment is counted once, from its forward-strand read.
  - Kept: proper pairs, mate mapped, edit distance ≤ 10.
  - **No MAPQ filter and no duplicate removal.**
- **Fragment centre:** even fragment lengths go to the left centre base in `bincount`; `NCPpeak` splits them 0.5/0.5.
- **Peak caller:** greedy and unsmoothed. Tied dyad counts are ordered exactly as Python 2 sorted them. The original's comparator never returns "equal", and Python 3.6 and 3.12 were checked to reproduce Python 2.7's order through `cmp_to_key`. `--fast` uses a plain tie rule instead.
- **Text output:** numbers are rounded and printed as Python 2 did (round half away from zero; `%.12g`).
- **Not ported:** `bowtie-submit_ver2`, which the GEO driver calls, is missing from Park's repo. The align step uses `bowtie-submit`'s command with the driver's `trim_to = 50`.

## Departures (opt-in, config keys)

| Key | Default (= Park) | Why you might change it |
|---|---|---|
| `reads.min_mapq` | 0 (no filter) | multi-mappers are placed at random; e.g. 10 or 30 removes them |
| `reads.dedup` | false | drop fragments with identical start and length |
| `score.missing_as_nan` | false | write bins without input reads as NA instead of dropping them |
| `peaks --fast` (CLI only) | off | faster tie order; can pick different dyads among equal counts |

Sliding-window (`Bdata`) areas are computed, but not scored by default (GEO didn't publish them).
Run `python -m condseq score` on them if needed.

## Validation (synthetic data, 2026-10-07)

`tests/simulate.py` builds a 600 kb genome, as shown in this tree:
```
synthetic genome (600 kb)
├── 3 chromosomes
├── a duplicated 2 kb segment     (multi-mappers)
└── a 5 kb N gap                  (bins without reads)
3,127 nucleosomes, each with a known c½ in 60 kb domains
reads: paired 50 bp
├── shallow series: 6 titration points, 40k pairs each
└── deep: 3 points, 150k pairs each
```

**1. Against Park's original scripts** (`tests/equivalence/`, job 51109796). Park's Python 2
scripts were run on the same BAMs with the GEO run's arguments:

| Stage | Files | Rows identical | Text-identical cells |
|---|---|---|---|
| bincount | 18 | 18/18 | 17,820/17,820 |
| score (bins) | 14 | 14/14 | 4,585/4,585 |
| coverage | 9 | 9/9 | 1,758,580/1,758,580 |
| ncov | 9 | 9/9 | 9,381/9,381 |
| bindata | 9 | 9/9 | 71,307/71,307 |
| score (nucleosomes) | 6 | 6/6 | 6,254/6,254 |
| peaks | 3 | 3/3 | 3,127/3,127 |
| num | 6 | 6/6 | 360/360 |
| fit (4PL) | 1 | 1/1 | 420/420 |

**2. Against the simulated truth** (`tests/check_truth.py`):

| Output | Pearson r with truth |
|---|---|
| 10 kb scores | 0.98–0.998 (points 2–5); 0.80 at point 1, where only 8% of nucleosomes have left |
| 1 kb deep scores | 0.96–0.99 |
| Per-nucleosome scores | 0.89–0.95; 3,118 of 3,127 called dyads within 20 bp of a true one |
| c½ (10 kb) | 0.991, mean error 0.011 mM |

**Not validated yet:** real reads, because we have none. The plan for the first real data is
below, under "Validate on GEO".

Re-run the validation:
```bash
cd ~/condense-seq/condense-seq-pipeline
sbatch tests/run_synthetic.sbatch                   # ~1 min of work
sbatch tests/equivalence/run_original.sbatch        # ~2 min; prints the comparison table
mamba activate condseq && python -m pytest tests/test_units.py
python tests/check_truth.py ../lab/outputs/condense-seq-pipeline/synthetic
```

## Running on real data

1. **Env** (once): `sbatch slurm/00_make_env.sbatch` builds `condseq` on holylabs. It's already built.
2. **References** (once): `sbatch slurm/10_reference.sbatch` puts hg38 and mm10 with bowtie2 indexes
   in `../lab/data/ref/` (~8 GB). **Not run yet**, because lab storage is nearly full.
   Park used the 4DN references; see the script's header.
3. **Sample sheet** (`samples.tsv`, tab-separated), one row per library:
   ```
   cell  sample  agent  tnum  rep  depth    fastq1                      fastq2
   H1    NCP     sp     0     1    deep     fastq/a_R1.fastq.gz         fastq/a_R2.fastq.gz
   H1    NCP     sp     8     1    deep     fastq/b_R1.fastq.gz         fastq/b_R2.fastq.gz
   H1    NCP     sp     1     1    shallow  fastq/c_L1_R1.fastq.gz,fastq/c_L2_R1.fastq.gz  ...
   ```
   - `depth` is `deep` or `shallow`.
   - Each series (cell, sample, agent, rep, depth) needs its input, `tnum 0`.
   - List several lanes comma-separated.
4. **Titration table** per (cell, sample, agent): `titration/<cell>_<sample>_<agent>.tsv` with columns
   `tnum  conc  survival_fraction` (tnum 0 → conc 0, fraction 1).
   - The survival fraction is the measured soluble fraction (UV-Vis) at that point.
   - **It sets the score offset (−ln f), so it must come from the experiment.**
   - Park's 10-column format is also read.
   - **For the GEO samples,** Park's titration measurements are in the paper's Source Data (Fig. 1b,
     ED Figs. 3d, 7c, 8a, 10a). `tools/titration_from_source_data.py` turns them into tables for all
     20 published series, already in `../lab/data/park2025/titration/`.
5. **Config**: copy `config/config.example.yaml` next to the sample sheet, set `genomes`, `cells`, paths.
6. **Run** (from `condense-seq-pipeline/`, in the `condseq` env):
   ```bash
   snakemake -s workflow/Snakefile --configfile /path/config.yaml \
             --directory ../lab/outputs/condense-seq-pipeline/<run> --profile workflow/profiles/slurm -n   # dry run
   # then the same without -n; Snakemake resumes where it stopped
   ```
   Rules ask for 8–24 GB and up to 24 h (alignment). The heaviest step on deep human data is
   `peaks`, about 20 M candidate dyads on chr1.

### Validate on GEO (first thing to do with real reads)

`geo_test/` re-runs GEO samples from their raw reads and compares the tables. The manifest,
download script, sample sheets, configs and comparison script are there, in three steps:
- 1.3 GB: H1 PEG, 10 kb
- 8.8 GB: H1 spermine, shallow, 10 kb + c½
- ~57 GB: H1 spermine, deep, per nucleosome

Prepared 2026-10-12; nothing downloaded yet. See `geo_test/README.md`.

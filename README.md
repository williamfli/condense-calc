# condseq

Turns condense-seq reads (paired-end FASTQ) into condensability tables and a QC report.
It reimplements the processing behind the Park et al. 2025 GEO tables (GSE252941); on test data its
output is identical to the authors' scripts. The output uses the same table format.

## Install

```bash
git clone git@github.com:williamfli/condense-seq-pipeline.git
cd condense-seq-pipeline
mamba env create -f environment.yml      # bowtie2, samtools, snakemake, python packages
mamba activate condseq
pip install -e .
```

## Run

```bash
condseq reference hg38 ~/refs            # once: hg38 FASTA + bowtie2 index (~7 GB)
condseq init myexp                       # writes myexp/config.yaml, samples.tsv, titration/
# edit the three files (below), then:
condseq run myexp/config.yaml --cores 8  # results in myexp/results/
```

On a Slurm cluster, each step can run as its own job. The controller is light but runs for the
whole workflow, so start it as a small job too:

```bash
sbatch -c 1 --mem 4G -t 3-00:00 --wrap \
  "condseq run myexp/config.yaml --slurm --account mylab --partition shared"
```

`condseq run ... -n` shows what would run. Re-running the same command resumes after an
interruption.

## Inputs

**`samples.tsv`**: one row per library, tab-separated.

| column | meaning |
|---|---|
| `cell`, `sample`, `agent` | e.g. `H1`, `NCP`, `sp` (spermine); `cell` must be listed under `cells:` in the config |
| `tnum` | titration point; `0` is the input |
| `rep` | replicate number |
| `depth` | `deep` (also gives per-nucleosome scores) or `shallow` |
| `fastq1`, `fastq2` | read 1 and read 2; several lanes comma-separated |

A *series* is one cell, sample, agent, rep and depth. Each series needs its input (`tnum` 0).

**`titration/<cell>_<sample>_<agent>.tsv`**: the measured soluble fraction at each point.

```
tnum  conc   survival_fraction
0     0      1.0
8     0.794  0.105
```

`survival_fraction` (*f*) is the fraction of nucleosomes left in the supernatant (e.g. by UV-Vis).
It sets the score offset (−ln *f*), so it must come from the experiment. `conc` is used for the
c½ fit.

**`config.yaml`**: where the samples, titration tables and genomes are, and which genome each cell
uses. Relative paths are read relative to the config file.

## Outputs (`results/`)

| file | content |
|---|---|
| `tables/<series>_1kb_score.gtab.gz`, `_10kb_` | condensability score per bin, one column per titration point |
| `tables/<series>_score/<series>_<chr>_score.gtab.gz` | per-nucleosome scores (deep series) |
| `tables/<series>_10kb_num.gtab.gz`, `_10kb_4PL.gtab.gz` | nucleosome numbers and c½ fit (series with ≥ 5 points) |
| `qc/report.html` | QC figures and tables; `qc/flags.tsv` lists anything outside the expected range |
| `bam/` | alignments (can be deleted once the tables are made) |

Score at a titration point with soluble fraction *f*, for supernatant count *s* and input count *c*:
−ln[(s+1)/Σ(s+1) ÷ (c+1)/Σ(c+1)] − ln *f*. Higher means more condensable.

## QC

Set `qc:` in the config, or `--qc` on the command line:

| level | contents |
|---|---|
| `off` | none |
| `default` | reads, alignment rate, multi-mapping share, fragment lengths, titration, score distributions, score vs input and supernatant counts, c½ fit success, nucleosome count and spacing (deep), replicate agreement |
| `verbose` | adds duplicate rate, GC dependence, a chromosome-wide score profile, dinucleotide periodicity at nucleosome dyads |

Several checks follow the paper's Extended Data Fig. 2 (fragment lengths, input vs supernatant,
dinucleotide periodicity, replicate correlation). Flags are warnings only and never stop a run.
To change a threshold:

```yaml
qc_thresholds: {min_replicate_rho: 0.7, max_lowmapq_fraction: 0.2}
```

The thresholds and their meanings are listed at the end of the report.

## Settings

The defaults reproduce the GEO tables. These keys can be added to the config:

| key | default | |
|---|---|---|
| `threads` | 8 | bowtie2 threads per sample |
| `bin_sizes` | `[1000, 10000]` | |
| `reads: {min_mapq: 0}` | no MAPQ filter, as Park | e.g. 10 drops multi-mapping reads, whose placement depends on the genome index |
| `reads: {dedup: false}` | duplicates kept, as Park | |
| `reads: {mm: 10}` | max. edit distance | |
| `score: {missing_as_nan: false}` | bins without input reads are dropped | `true` writes them as NA |
| `mscale` | `{Ca: 5}`, others 1 | molecule-number scale per agent |
| `fit` | sigmoid, ≥ 5 points, 10 kb, R² ≥ 0.5 | |
| `trim_to` | 50 | reads are trimmed to 50 bp before alignment |

## Tests

```bash
pytest tests/test_units.py
bash tests/synthetic/run.sh /tmp/condseq-test      # simulated experiment, ~1 min
python tests/synthetic/check_truth.py /tmp/condseq-test
```

`tests/equivalence/run.sh` runs the authors' original Python 2 scripts on the same simulated reads
and compares every table. See [docs/details.md](docs/details.md) for the processing steps and
validation.

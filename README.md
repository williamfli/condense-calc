# condense-calc

Turns condense-seq reads (paired-end FASTQ) into condensability scores and a QC report.
The processing is that of Park et al. 2025 (GEO GSE252941); on test data the output is identical
to the authors' scripts.

## Install

```bash
git clone git@github.com:williamfli/condense-calc.git
conda env create -f condense-calc/environment.yml     # installs everything, ~10 min
```

## Run

```bash
conda activate condense-calc
condense init path/to/fastqs        # writes samples.tsv, one row per FASTQ
# fill in the conc and fraction columns, then:
condense samples.tsv                # results in results/ next to samples.tsv
```

On a Slurm cluster, add `--slurm` (and `--account`, `--partition` if needed). `condense` then
submits itself as a small job, which runs each step as its own job, and returns; this also works
from an interactive session (salloc, Open OnDemand), which can then be closed. `--mail ADDRESS`
sends an email if the run fails. Without `--slurm`, `condense` runs on the machine it is started on,
with the CPUs and memory it may use there (in a Slurm job: that job's allocation); on a cluster login
node it stops and asks for `--slurm` (or `--cores N` to run there anyway).

The first run downloads the genome (hg38 by default, ~7 GB, ~30 min) into `~/.cache/condense-calc`;
with `--slurm` this, or building the index of `--genome FASTA`, runs as its own job first.
To keep genomes elsewhere or share them, set `CONDENSE_REFS=/path/to/refs`.

Re-running the same command resumes an interrupted run. `-n` checks the sample sheet and shows what
would run (it does not download a genome).

## Sample sheet

```
series   conc    fraction  fastq
H1_sp    0       1         reads/input_R1.fastq.gz
H1_sp    0.794   0.105     reads/sp8_L001_R1.fastq.gz
H1_sp    0.794   0.105     reads/sp8_L002_R1.fastq.gz
```

| column | |
|---|---|
| `series` | one titration; each series needs an input row (`conc` 0) |
| `conc` | agent concentration (any unit; used for the c½ fit) |
| `fraction` | soluble fraction measured at that concentration (input: 1). It sets the score offset, −ln *fraction* |
| `fastq` | read-1 file; read 2 is found by name (`_R1` → `_R2`, `_1.fastq` → `_2.fastq`). Wildcards work. Rows with the same series and conc are lanes of one library |
| `rep` | optional: replicate number; replicates of a series are compared in the QC report |
| `fastq2`, `agent` | optional: read-2 file if its name doesn't follow read 1; agent name (`Ca` scales nucleosome numbers by 5) |

## Results

| file | |
|---|---|
| `tables/<series>_1kb_score.gtab.gz`, `_10kb_` | condensability score per bin, one column per concentration |
| `tables/<series>_10kb_4PL.gtab.gz` | c½ per bin (series with ≥ 5 concentrations) |
| `tables/<series>_score/` | per-nucleosome scores (with `--nucleosomes`) |
| `qc/report.html` | QC figures; `qc/flags.tsv` lists values outside the expected range |
| `bam/` | alignments; can be deleted afterwards |

Tables are gzipped tab-separated text, in the GEO format. Score for supernatant count *s* and input
count *c* in a bin: −ln[(s+1)/Σ(s+1) ÷ (c+1)/Σ(c+1)] − ln *fraction*. Higher is more condensable.

## Options

| option | |
|---|---|
| `-g, --genome` | `hg38` (default), `mm10`, or a FASTA file (indexed on first use) |
| `--nucleosomes` | also score every nucleosome; for deeply sequenced series |
| `--qc off\|default\|verbose` | QC report detail (verbose adds duplicates, GC, dinucleotide periodicity) |
| `--slurm`, `--account`, `--partition` | run on Slurm |
| `--mail ADDRESS` | Slurm: email if the run fails or reaches 80% of its time limit |
| `-o` | output folder |
| `--settings FILE` | YAML with the settings below |

Settings (defaults reproduce the GEO tables):

| key | default | |
|---|---|---|
| `reads: {min_mapq: 0}` | no filter | e.g. 10 drops multi-mapping reads |
| `reads: {dedup: false}` | duplicates kept | |
| `bin_sizes` | `[1000, 10000]` | |
| `threads` | 8 | bowtie2 threads per library |
| `qc_thresholds` | see `qc/report.html` | e.g. `{min_replicate_rho: 0.7}` |

More settings, multi-genome runs and the processing steps: [docs/details.md](docs/details.md).

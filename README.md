# condense-calc

condense-calc turns condense-seq reads into condensability scores. You give it paired-end FASTQs of the input and the supernatants of a titration, with each concentration and its soluble fraction; it gives scores per 1 kb bin, 10 kb bin and nucleosome, c½ per 10 kb bin, and a QC report.
It is based on the processing of Park & Ha, Nature, 2025, with a few design changes; `--park` reproduces their processing ([docs/park.md](docs/park.md)).

## Install

```bash
git clone git@github.com:williamfli/condense-calc.git
conda env create -f condense-calc/environment.yml
conda activate condense-calc
```

## Your data

Put both reads of every library in one folder. Read 2 has the name of read 1 with `_R2` for `_R1` (or `_2` for `_1`).

```
my_experiment/fastq/
├── input_R1.fastq.gz      input_R2.fastq.gz
├── sp8_L001_R1.fastq.gz   sp8_L001_R2.fastq.gz
└── sp8_L002_R1.fastq.gz   sp8_L002_R2.fastq.gz
```

```bash
cd my_experiment
condense init fastq        # writes samples.tsv, one row per FASTQ: fill in series, conc and fraction
```

```
series   conc    fraction  fastq
H1_sp    0       1         fastq/input_R1.fastq.gz
H1_sp    0.794   0.105     fastq/sp8_L001_R1.fastq.gz
H1_sp    0.794   0.105     fastq/sp8_L002_R1.fastq.gz
```

- `series`: name of one titration (letters, digits, `_`, `-`). Each series needs an input row.
- `conc`: concentration of the agent, any unit. 0 for the input.
- `fraction`: soluble fraction measured at that concentration, not rounded. 1 for the input.
- `fastq`: read-1 file. Rows with the same series and conc are lanes of one library.
- Optional columns: `rep` (replicate number), `agent` (`Ca` scales nucleosome numbers by 5), `fastq2` (read-2 file).

## Run

```bash
condense samples.tsv                                     # on this computer
condense samples.tsv --slurm --account ACCOUNT --partition PARTITION --mail ADDRESS   # on a Slurm cluster
condense samples.tsv -g mm10                             # mouse (default: -g hg38)
condense samples.tsv -g genome.fa                        # another genome (FASTA file)
condense samples.tsv --nucleosomes                       # a sheet of deep series: also a score per nucleosome
```

Rerun the same command to resume. The first run downloads the genome to `~/.cache/condense-calc` (set `CONDENSE_REFS` to use another folder).

## Results

| file, in `results/` next to `samples.tsv` (or `-o FOLDER`) | |
|---|---|
| `tables/<series>_10kb_score.gtab.gz`, `_1kb_` | score per bin, one column per concentration; higher is more condensable |
| `tables/<series>_10kb_4PL.gtab.gz` | c½ per bin (series with 5 or more concentrations, input included) |
| `tables/<series>_score/` | score per nucleosome, one file per chromosome (`--nucleosomes`) |
| `qc/report.html` | QC report; `qc/flags.tsv` lists values outside the expected range |
| `run_info.json` | settings, references and software versions of the run |

## Reproduce Park et al.

```bash
condense fetch SRR31626986 SRR31626987 SRR31626988 SRR31626989 -o fastq   # e.g. H1 PEG rep 2, from NCBI SRA
condense init fastq                                 # then fill in samples.tsv as above
condense samples.tsv --park
```

What it takes to get the GEO tables byte for byte: [docs/park.md](docs/park.md).

Settings, config files, steps, resources and tests: [docs/details.md](docs/details.md).

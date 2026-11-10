# condense-calc

condense-calc turns condense-seq reads into condensability scores.
- **You give:** paired-end FASTQs of the input and of the supernatant at each concentration of a titration, and the soluble fraction at each concentration.
- **You get**, for every series: a condensability score per 10 kb and per 1 kb bin at each concentration, a QC report, and, for series with the input and at least 4 concentrations, c½ per 10 kb bin: the concentration at which half of the bin's nucleosomes stay soluble.
- **Deep series** (~300 M read pairs per library) can also get a score per nucleosome. Add `--nucleosomes` to the run command.

By default it uses the processing of the condense-seq paper (Park, Merino-Urteaga, Karwacki-Neisius, Carrizo, Athreya … Ha, Nature, 2025): from the paper's reads, the tables are identical to the paper's GEO tables. `--revised` makes a few design changes ([what and why](docs/revised.md)).

## Install

```bash
git clone https://github.com/williamfli/condense-calc.git
cd condense-calc
conda env create -f environment.yml
conda activate condense-calc
```

To update later, in `condense-calc`: `git pull && conda env update -f environment.yml`. Linux and macOS; on a Mac with Apple silicon the default mode's bowtie2 2.4.1 runs under Rosetta 2 (if it is missing, condense-calc says how to install it).

## Your data

Each experiment gets its own folder, separate from the `condense-calc` folder you installed from. In it, a `fastq` folder holds two files per library (paired-end): read 1 (`..._R1.fastq.gz`) and read 2 (`..._R2.fastq.gz`, the same name with `_R2` for `_R1`, or `_2` for `_1`). Results are written next to them:

```
my_experiment/
├── fastq/
│   ├── input_R1.fastq.gz      input_R2.fastq.gz
│   ├── sp8_L001_R1.fastq.gz   sp8_L001_R2.fastq.gz
│   └── sp8_L002_R1.fastq.gz   sp8_L002_R2.fastq.gz
├── samples.tsv                written by condense init; you fill it in
└── results/                   written by condense
```

1. Make the folders and put the FASTQ files in `fastq` (copy, move or download them there):

   ```bash
   mkdir -p ~/my_experiment/fastq
   cp /path/to/your/reads/*.fastq.gz ~/my_experiment/fastq/
   ```

2. In a new terminal, activate the environment (once per terminal) and go into the experiment folder:

   ```bash
   conda activate condense-calc
   cd ~/my_experiment
   ```

3. Write the sample sheet, one row per read-1 file. `--series` names the titration:

   ```bash
   condense init fastq --series H1_sp
   ```

4. Open `samples.tsv` in a text editor or a spreadsheet and fill in `conc` and `fraction` (keep the tabs between columns; from a spreadsheet, save as tab-separated text):

   ```
   series   conc    fraction  fastq
   H1_sp    0       1         fastq/input_R1.fastq.gz
   H1_sp    0.794   0.105     fastq/sp8_L001_R1.fastq.gz
   H1_sp    0.794   0.105     fastq/sp8_L002_R1.fastq.gz
   ```

   - `series`: name of one titration (letters, digits, `_`, `-`). Each series needs an input row.
   - `conc`: concentration of the agent, any unit. 0 for the input.
   - `fraction`: soluble fraction measured at that concentration, with all the digits you have. 1 (or empty) for the input.
   - `fastq`: read-1 file. Rows with the same series, rep and conc are lanes of one library.
   - Optional columns: `rep` (replicate number; the series is then named `<series>_rep<N>`), `agent` (e.g. `PEG`, `Ca`), `fastq2` (read-2 file, if not named as above).

## Run

From the experiment folder, with the environment active:

```bash
condense samples.tsv -n                 # check the sheet and show what would run
condense samples.tsv                    # run on this computer
condense samples.tsv --slurm            # run on a Slurm cluster (optional: --account A --partition P --mail ADDRESS)
condense samples.tsv -g mm10            # mouse (default: -g hg38)
condense samples.tsv --revised          # with the design changes of docs/revised.md
condense samples.tsv -g genome.fa       # another genome (FASTA file)
condense deep.tsv --nucleosomes         # a sheet of deep series only (~300 M read pairs per library): also a score per nucleosome
condense deep.tsv --nucleosomes --states auto   # also QC per chromatin state (built in: GM12878, E14 mESC; else a BED file)
```

Neither the species nor the sequencing depth is detected from the reads. Give `-g` for any genome other than human, and `--nucleosomes` for deep series. Both apply to every series in the sheet, so put mouse and human series, or deep and shallow ones, in separate sheets.

Rerun the same command to resume. The first run downloads the genome to `~/.cache/condense-calc` (set `CONDENSE_REFS` to use another folder) and installs bowtie2 2.4.1 there. For hg38 it also builds the bowtie2 index (~1.5 h with 8 cores; with `--revised`, NCBI's prebuilt index, ~7 GB). With `-g genome.fa` it builds an index next to the FASTA (hours for a mammalian genome).

## Results

In `results/` next to `samples.tsv` (or in `-o FOLDER`). Tables are gzipped tab-separated text.

| file | contents |
|---|---|
| `tables/<series>_10kb_score.gtab.gz`, `tables/<series>_1kb_score.gtab.gz` | score per bin; higher is more condensable. One column per concentration, input excluded: `<series>_1`, `<series>_2`, ... in increasing order (as in `titration/<series>.tsv`). Bins without input reads are left out (`--revised`: NA) |
| `tables/<series>_10kb_4PL.gtab.gz` | c½ per bin: the concentration at which half of the bin's nucleosomes stay soluble (series with the input and at least 4 concentrations) |
| `tables/<series>_score/` | score per nucleosome, one file per chromosome (`--nucleosomes`) |
| `qc/report.html` | QC report ([example](docs/qc_example.md)), with panels drawn as in the paper's Extended Data Fig. 2; `qc/flags.tsv` lists values outside the expected range |
| `run_info.json` | settings, references, condense-calc and bowtie2 versions |

## Reproduce the paper's GEO tables

```bash
condense fetch SRR31626986 SRR31626987 SRR31626988 SRR31626989 -o fastq   # H1 PEG rep 2, from NCBI SRA
condense init fastq --series H1_NCP_PEG
```

Fill in `samples.tsv` (conc and fraction from the paper's Source Data):

```
series      conc  fraction     fastq
H1_NCP_PEG  9.5   0.097744361  fastq/SRR31626986_1.fastq.gz
H1_NCP_PEG  9.5   0.097744361  fastq/SRR31626987_1.fastq.gz
H1_NCP_PEG  0     1            fastq/SRR31626988_1.fastq.gz
H1_NCP_PEG  0     1            fastq/SRR31626989_1.fastq.gz
```

```bash
condense samples.tsv
```

On a cluster, run `condense fetch` in a job (it prints an `sbatch` line to paste) and add `--slurm` to the last line.
Name GM12878 and mCD8T series as on GEO (`GM_...`, `mCD8T_...`), so chrY is left out as in the paper.
For GEO's file and column names too: [docs/geo.md](docs/geo.md).

Settings, config files, steps, resources and tests: [docs/details.md](docs/details.md).

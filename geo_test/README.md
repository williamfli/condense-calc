# First real-data check: re-run GEO samples and compare with the GEO tables

**Status (2026-10-12): prepared, nothing downloaded yet.** Bill approved steps 0 and 1 with raw
data on netscratch. The downloads and the hg38 reference wait for his go-ahead.

| Step | Series (ENA runs in `manifest.tsv`) | Download | Read pairs | Compare with |
|---|---|---|---|---|
| 0 | H1 NCP PEG rep2: input + ID 6, 2 runs each | 1.28 GB | 21.3 M | `GSE252941_H1_NCP_PEG_2rep_10kb_score` |
| 1 | H1 NCP spermine rep2, shallow: input + IDs 1–9, 2 runs each | 8.77 GB | 146.9 M | `GSE252941_H1_NCP_sp_2rep_10kb_score`; c½ fit |
| 2 (later) | H1 NCP spermine rep1, deep: input, IDs 4 and 8 | ~57 GB | ~1 B | deep 1 kb and per-nucleosome tables |

The runs come from the ENA file report for PRJNA1063283, fetched 2026-10-12. The manifest
records each run's URL, md5 and size.

## Inputs
- **FASTQ:** `/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/raw/GSE252941/` (`download.sbatch`, md5-checked).
  netscratch deletes files after 90 days; download again if needed.
- **Titration:** `../../lab/data/park2025/titration/H1_NCP_{PEG,sp}.tsv`. These are made from Park's
  Source Data (Fig. 1b, ED Fig. 8a) by `tools/titration_from_source_data.py`; they are Park's
  own titration measurements.
  - PEG ID 6: 9.5% PEG, f = 0.098.
  - Spermine ID 8: 0.794 mM, f = 0.105.
- **Reference:** hg38 from `../slurm/10_reference.sbatch` (~8 GB in `../../lab/data/ref/`). Park used
  the 4DN hg38 and a bowtie2 version we don't know, so exact agreement isn't expected.

## Running
```bash
cd ~/condense-seq/condense-seq-pipeline
sbatch slurm/10_reference.sbatch                          # once, ~8 GB
cd geo_test && sbatch download.sbatch 0 && cd ..          # or: download.sbatch 0 1
module load python/3.12.11-fasrc02 && mamba activate condseq
RUN=/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/pipeline_runs/geo_step0
snakemake -s workflow/Snakefile --configfile geo_test/config_step0.yaml --directory $RUN --profile workflow/profiles/slurm
python geo_test/compare_geo.py $RUN/tables/H1_NCP_PEG_2rep_10kb_score.gtab.gz \
       ../lab/data/GSE252941/GSE252941_H1_NCP_PEG_2rep_10kb_score.gtab.txt.gz
```
BAMs and intermediate files stay in the netscratch run folder. Copy `tables/` to
`../lab/outputs/condense-seq-pipeline/geo_step<N>/` to keep the results.

## What would count as success
Rows should be nearly the same set, apart from bins at assembly gaps and on chrY, where the
reference build matters. We'd expect r ≥ 0.99 and a mean difference near 0 for every column.
A constant offset would point to the titration value. A spread that's otherwise high-correlation
would point to the reference or aligner.

# First real-data check: re-run GEO samples and compare with the GEO tables

**Status (2026-10-13): steps 0 and 1 done (step 0 write-up: `step0/README.md`).** hg38 is installed
in `../../lab/data/ref/hg38` (NCBI no-alt analysis set + its bowtie2 index, ~7 GB).

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

## Step 0 result (2026-10-13)

Full write-up, with inputs, every step, figures and the scripts: **[`step0/README.md`](step0/README.md)**.

- **No offset:** mean difference −0.0003, so the titration value (f = 0.098) and the normalization match.
- **Identical counts:** 63% of bins have the same fragment counts as Park's.
- **Agreement:** r = 0.964 over all bins. In bins with no MAPQ < 10 fragments (58% of bins),
  r = 0.995, and 87% of those bins have identical counts.
- **Where the differences are:** they follow the share of multi-mapping fragments in a bin and sit
  in centromeric and other repeats.
- **Not the aligner:** bowtie2 2.3.5.1, 2.4.5 and 2.5.5, and seeds 0 and 1, all give identical
  counts. A multi-mapper's placement is a fixed function of the read and the index.
  - The likely cause is Park's different index (4DN "4D_hg38"). 4DN needs an account key to
    download it, so this is untested.
- **Scale:** the ours-vs-GEO differences (sd 0.10) are small next to the difference between Park's
  two replicates (sd 0.41).

Run folder: `/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/pipeline_runs/geo_step0` (BAMs,
intermediates, realignments). Tables kept in `../../lab/outputs/condense-seq-pipeline/geo_step0/`.

## Step 1 result (2026-10-13)

H1 NCP spermine, replicate 2, shallow: input + IDs 1–9, 146.9 M read pairs. Same settings as step 0.

| Column | Bins in both | Pearson | Spearman | Mean diff | Median abs. diff | abs. diff < 0.05 |
|---|---|---|---|---|---|---|
| `sp_1` (0.158 mM) | 292,549 | 0.940 | 0.955 | 0.0000 | 0.0006 | 88.4% |
| `sp_2` | 292,549 | 0.946 | 0.961 | 0.0001 | 0.0012 | 88.0% |
| `sp_3` | 292,549 | 0.956 | 0.973 | 0.0000 | 0.0011 | 88.2% |
| `sp_4` (0.316 mM) | 292,549 | 0.962 | 0.979 | −0.0000 | 0.0011 | 87.9% |
| `sp_5` | 292,549 | 0.973 | 0.985 | 0.0002 | 0.0011 | 88.1% |
| `sp_6` | 292,549 | 0.974 | 0.986 | 0.0001 | 0.0011 | 87.7% |
| `sp_7` | 292,549 | 0.974 | 0.986 | 0.0003 | 0.0002 | 87.4% |
| `sp_8` (0.794 mM) | 292,549 | 0.973 | 0.985 | 0.0005 | 0.0009 | 86.3% |
| `sp_9` (1 mM) | 292,549 | 0.977 | 0.987 | 0.0003 | 0.0012 | 87.6% |

- **No offset at any point,** so Park's soluble fraction is right for all nine titration points.
- **The same pattern as step 0:** 86–88% of bins within 0.05 of GEO, and the same bins missing
  (1,402 only in GEO, 16 only in ours). r is lowest at the first points, where scores vary least
  between bins.
- **c½ fit (10 kb):** 292,373 of 292,565 bins fit with R² ≥ 0.5 (median R² 0.93). Median c½ is
  0.215 mM, in line with the genome-wide titration (0.48 surviving at 0.200 mM). GEO has no c½
  table, so this is a sanity check, not a comparison.
- **1 kb table:** built, but GEO has no shallow 1 kb spermine table to compare it with.
- **Download:** ENA served 14 of 40 files at full length with a few hundred bytes wrong;
  `download.sbatch` now downloads a file again until its md5 matches.

Run folder: `/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/pipeline_runs/geo_step1`. Tables kept
in `../../lab/outputs/condense-seq-pipeline/geo_step1/`. Compared with
`python geo_test/compare_geo.py <ours> ../lab/data/GSE252941/GSE252941_H1_NCP_sp_2rep_10kb_score.gtab.txt.gz`.

## What would count as success
Rows should be nearly the same set, apart from bins at assembly gaps and on chrY, where the
reference build matters. We'd expect r ≥ 0.99 and a mean difference near 0 for every column (step 0: met in uniquely
mapping bins only; see above).
A constant offset would point to the titration value. A spread that's otherwise high-correlation
would point to the reference or aligner.

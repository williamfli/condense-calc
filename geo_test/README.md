# First real-data check: re-run GEO samples and compare with the GEO tables

**Status (2026-10-12): step 0 done (results below); step 1 not downloaded yet.** hg38 is installed
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

## Step 0 result (2026-10-12)

Run: `/n/netscratch/zhuang_lab/Lab/wfli/condense-seq/pipeline_runs/geo_step0` (BAMs, intermediates);
tables and the MAPQ counts kept in `../../lab/outputs/condense-seq-pipeline/geo_step0/`.
Input 12.7 M fragments, ID 6 7.9 M (after Park's read filters); alignment took ~15 min per sample.
GEO has no 1 kb PEG table, so only the 10 kb table is compared.

| | Bins | Pearson r | |d| < 0.05 | |d| > 0.5 |
|---|---|---|---|---|
| All bins in both tables | 292,544 | 0.964 | 83% | 0.9% |
| Bins with no low-MAPQ (< 10) fragments | 171,058 (58%) | 0.995 | 97% | 0.1% |
| Bins with > 50% low-MAPQ fragments | 15,184 (5%) | | median |d| 0.19 | 14% |

- Mean difference −0.0003, median |d| 0.0005: no offset, so the titration value (f = 0.098) and
  the normalization match.
- Differences follow the share of multi-mapping fragments in a bin (`mapq_bins.py`), not depth.
  The worst bins are pericentromeric and segmental-duplication regions (chr15 18–21 Mb,
  chr20 28–30 Mb, chrX 58–59 Mb). Park kept all alignments regardless of MAPQ, so bowtie2 places
  multi-mappers at random among equal hits; a different reference copy or bowtie2 version
  places them differently. This can't be removed while matching Park's settings.
- At ~25 supernatant fragments per 10 kb bin, a few moved fragments change a score by 0.1–0.5.
- Rows: 13 only ours, 1,408 only GEO (bins whose input count differs between 0 and > 0).

Conclusion: the pipeline reproduces GEO where reads map uniquely; elsewhere it is limited by
random multi-mapper placement. `reads.min_mapq` (a departure from Park) removes those reads.

## What would count as success
Rows should be nearly the same set, apart from bins at assembly gaps and on chrY, where the
reference build matters. We'd expect r ≥ 0.99 and a mean difference near 0 for every column (step 0: met in uniquely
mapping bins only; see above).
A constant offset would point to the titration value. A spread that's otherwise high-correlation
would point to the reference or aligner.

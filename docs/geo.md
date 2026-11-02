# Reproducing our GEO tables

The default processing is that of our paper (Park, Merino-Urteaga, Karwacki-Neisius, Carrizo, Athreya … Ha, Nature, 2025). From our reads it gives tables identical to GSE252941's.
Four things have to match:

1. **Reads as we aligned them.** `condense fetch` downloads runs from NCBI SRA with full base
   qualities and the instrument's read names. ENA's FASTQs differ: ENA serves 92 of the 406 GSE252941 runs as
   SRA Lite (all qualities `?`), and its read names change where bowtie2 puts multi-mapping reads.
2. **The default processing:** UCSC hg38.fa, bowtie2 2.4.1, MAPQ 0, bins without input reads dropped,
   Python 2's tie order ([revised.md](revised.md) lists these). Do not add `--revised`. A config file names its
   genome itself: point it at `hg38_ucsc` (the run warns if it is not the paper's reference).
3. **Soluble fractions at full precision.** Our titration files are not published; the paper's
   Source Data have the fractions at full precision.
4. **Names as on GEO.** With a [config file](details.md#config-files) (cell, sample, agent, tnum, rep and depth,
   `deep` for the deep series, as on GEO), file and column names match GEO's. With a sample sheet the numbers
   match under other column names; name series as on GEO (`GM_NCP_sp`, `mCD8T_WT_NCP_sp`) or give `chroms`,
   so chrY is left out where the paper left it out.

```bash
condense fetch SRR31626986 SRR31626987 SRR31626988 SRR31626989 -o fastq --threads 4
condense reference hg38_ucsc --bowtie2 2.4.1 --threads 8    # once, ~1.5 h; prints the FASTA and index paths
condense config.yaml
```

On a cluster, run the first two in a job (`condense` prints an `sbatch` line to paste) and add `--slurm` to the third.
`config.yaml` for this example (H1 PEG rep 2) is in [details.md](details.md#config-files).

Tables tested identical to GEO's (decompressed): `H1_NCP_PEG_2rep` (10 kb), `H1_NCP_sp_2rep` (10 kb),
`H1_NCP_sp_1rep_deep` and `H1_NCP_sp_2rep_deep` (1 kb and per nucleosome), `GM_NCP_sp_1rep` and `GM_NCP_sp_2rep`
(10 kb), `GM_NCP_sp_1rep_deep` and `GM_NCP_sp_2rep_deep` (1 kb and per nucleosome), `E14_NCP_sp_1rep` and
`E14_NCP_sp_2rep` (10 kb; GEO names these columns by library, so after renaming its header line).

# Reproducing the paper's GEO tables

The default processing is that of the condense-seq paper (Park, Merino-Urteaga, Karwacki-Neisius, Carrizo, Athreya … Ha, Nature, 2025). From the paper's reads it gives tables identical to GSE252941's.
Four things have to match:

1. **Reads as the paper aligned them.** `condense fetch` downloads runs from NCBI SRA with full base
   qualities and the instrument's read names. ENA's FASTQs differ: ENA serves 92 of the 406 GSE252941 runs as
   SRA Lite (all qualities `?`), and its read names change where bowtie2 puts multi-mapping reads.
2. **The default processing:** UCSC hg38.fa, bowtie2 2.4.1, MAPQ 0, bins without input reads dropped,
   Python 2's tie order ([revised.md](revised.md) lists these). Do not add `--revised`. A config file names its
   genome itself: point it at `hg38_ucsc` (the run warns if it is not the paper's reference).
3. **Soluble fractions at full precision.** The paper's titration files are not published; the paper's
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

Not identical: `mCD8T_DFMO_NCP_sp_2rep` (10 kb), checked for point 8 from SRR31626776, SRR31626777, SRR31627081 and
SRR31627082. Against GEO's column, r = 0.99975: 95% of bins agree within 1e-5, about 5% differ by roughly one read,
and GEO has one bin more. The same run on macOS (Apple silicon, bowtie2 2.4.1 under Rosetta 2) and on Linux gave
byte-identical bin counts and 10 kb and 1 kb tables, so the difference comes from the reads: the SRA runs of this
mouse CD8 series seem to differ slightly from those the paper aligned.

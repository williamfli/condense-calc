# Reproducing the paper's GEO tables

The default processing is that of the condense-seq paper (Park, Merino-Urteaga, Karwacki-Neisius, Carrizo, Athreya … Ha, Nature, 2025). From the paper's reads it gives tables identical to GSE252941's.
Four things have to match:

1. **Reads as the paper aligned them.** `condense fetch` downloads runs from NCBI SRA with full base
   qualities and the instrument's read names. ENA's FASTQs differ: ENA serves 92 of the 406 GSE252941 runs as
   SRA Lite (all qualities `?`), and its read names change where bowtie2 puts multi-mapping reads.
2. **The default processing:** UCSC hg38.fa, bowtie2 2.4.1, MAPQ 0, bins without input reads dropped,
   Python 2's tie order ([revised.md](revised.md) lists these). Do not add `--revised`. A config file names its
   genome itself: point it at `hg38_ucsc` (the run warns if it is not the paper's reference). Mouse: the `mm10`
   preset (ENCODE) for E14, but UCSC's mm10.fa for the mCD8T series (see [below](#macos-and-the-mcd8t-series)).
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

Tables tested identical to GEO's (decompressed), all run on a Linux Slurm cluster: `H1_NCP_PEG_2rep` (10 kb),
`H1_NCP_sp_2rep` (10 kb), `H1_NCP_sp_1rep_deep` and `H1_NCP_sp_2rep_deep` (1 kb and per nucleosome),
`GM_NCP_sp_1rep` and `GM_NCP_sp_2rep` (10 kb), `GM_NCP_sp_1rep_deep` and `GM_NCP_sp_2rep_deep` (1 kb and per
nucleosome), `E14_NCP_sp_1rep` and `E14_NCP_sp_2rep` (10 kb; GEO names these columns by library, so after renaming
its header line), and point 8 of `mCD8T_DFMO_NCP_sp_2rep` (10 kb, with UCSC's mm10.fa; see below).

## macOS and the mCD8T series

We tested on macOS manually (by Bill). Results agree to correlation 0.99975. The test: point 8 of
`mCD8T_DFMO_NCP_sp_2rep` (input SRR31626776 and SRR31626777, point 8 SRR31627081 and SRR31627082; `-g mm10`),
against GEO's 10 kb table. The same run on the Slurm cluster gave byte-identical bin counts and tables, so macOS
(bowtie2 2.4.1 under Rosetta 2) changes nothing.

The difference from GEO is the reference. The paper aligned the mCD8T series to UCSC's mm10.fa, not to the ENCODE
mm10 of the `mm10` preset (which reproduces E14). The two have the same sequences except on chrY, where ENCODE
hard-masks 98,854 more bases (chrY:90,745,846-90,844,699, at the pseudoautosomal region), and ENCODE orders the
contigs differently (chrUn before chrX). bowtie2 places a multi-mapping read among the hits it samples from the
read's 22-base seeds, so changing the index moves some highly repetitive reads: about 0.2% of the multi-mapping
fragments, 92-94% of them with a seed that also occurs in that chrY stretch. About 5% of the 10 kb bins then
differ by a read or two, the totals change, and every score shifts by about 1e-5. ENCODE's mm10 with UCSC's chrY
still differs in 7,686 bins, so both the mask and the contig order matter.

With UCSC's mm10.fa the table is identical to GEO's (all 256,090 rows):

```bash
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz && gunzip mm10.fa.gz
condense samples.tsv -g mm10.fa        # first builds the bowtie2 index next to mm10.fa (~40 min)
```

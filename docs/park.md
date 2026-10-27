# Differences from Park et al. (`--park`)

condense-calc follows the processing of Park & Ha, Nature, 2025 (GEO GSE252941; scripts:
[spark159/condense-seq](https://github.com/spark159/condense-seq)). By default it changes a few things.
`--park` (or `park: true` in a settings or config file) does them as Park et al. did.

## What differs

| what | condense-calc default | `--park` (Park et al.) | why |
|---|---|---|---|
| Reads that fit several places in the genome equally well (multi-mapping; `reads: {min_mapq}`) | dropped: a fragment counts only if its forward-strand read has MAPQ ≥ 10. Some 10 kb bins (mostly centromeres and other repeats) then have no input reads and are NA | kept (0: no filter) | Where bowtie2 puts such a read depends on the read names and base qualities, so their counts differ between pipelines. |
| Bins without input reads (`score: {missing_as_nan}`) | written as NA rows | dropped | Every bin table has every bin, so bin tables from different runs line up row for row. |
| Human reference (`-g hg38`) | NCBI GRCh38 no-alt analysis set | UCSC hg38.fa (`hg38_ucsc`: 455 sequences, 261 of them alt haplotypes) | Alt haplotypes are extra copies of parts of the main chromosomes, so bowtie2 sees reads there as multi-mapping. |
| bowtie2 (`bowtie2`) | 2.5.5, the environment's | 2.4.1 (Park et al.'s Reporting Summary), installed once as a conda env in the refs folder | The environment's version needs no second install. On the H1 libraries tested the two gave the same counts; on E14 rep 1, 2.5.5 differed by one count in 10 bins. |
| FASTQ check (`fastq_check`) | `lite`: stop on SRA Lite reads (base qualities all `?`) | `park`: also stop on read names not in the form Park et al. aligned (`@<spot> 1:N:0:<index>`) | Read names change where bowtie2 puts reads with several equally good hits; almost all of these have MAPQ < 10 and are dropped by default. |
| Nucleosome calling: positions with equal counts | taken in order of position | in the order Python 2's sort left them | Park's order is a side effect of Python 2's sort, not a rule. |
| Fragment centre | from the fragment's two ends | Park's formula | Park's formula is one base off for some odd-length fragments and wraps a centre before position 0 into the last bin. bowtie2's end-to-end paired alignment produces neither case, so no output changes. |
| chrY, when `chroms` is not given | kept: all main chromosomes | dropped for GM (GM12878) and mCD8T, matched by the first part of the series or cell name (`GM_...`, `mCD8T_...`) | Park's driver did this by name. The default does not guess from names; give `chroms` to drop chrY. |

Mouse (`-g mm10`) is the ENCODE mm10 no-alt analysis set in both modes: equivalent to the index Park et al. used
(the E14 tables reproduce exactly).

The same in both modes:
- Soluble fractions are used at full precision, exactly as written in the sample sheet or titration file.
- Misspelled or invalid settings stop the run with a message.
- Each run writes `run_info.json`: mode, settings, references, and the bowtie2 version of each library.
- Sliding-window areas of deep series (Park's `Bdata` files) are made only with `bdata: true`. No GEO table uses them.
- Chromosomes of human and mouse genomes are in GEO order (1..N, X, Y), whatever the order in the FASTA.

Kept exactly as Park et al.: the steps and their parameters ([details.md](details.md#steps)), the score formula
(with 1 added to every count), the c½ fit, no duplicate removal, at most 10 mismatched or gapped bases per read,
and GEO's file and column names (with a [config file](details.md#config-files)).

## Reproducing the GEO tables exactly

1. **Reads as Park et al. aligned them.** `condense fetch` downloads runs from NCBI SRA with full base
   qualities and the instrument's read names. ENA's FASTQs differ: ENA serves 92 of the 406 GSE252941 runs as
   SRA Lite (all qualities `?`), and its read names change where bowtie2 puts multi-mapping reads.
2. **`--park`:** UCSC hg38.fa, bowtie2 2.4.1, MAPQ 0, bins without input reads dropped, Python 2's tie order.
   A config file names its genome itself: point it at `hg38_ucsc` (Park mode warns if it is not Park's reference).
3. **Soluble fractions at full precision.** Park et al.'s titration files are not published; the paper's
   Source Data have the fractions at full precision.
4. **Names as on GEO.** With a [config file](details.md#config-files) (cell, sample, agent, tnum, rep and depth,
   `deep` for the deep series, as on GEO), file and column names match GEO's. With a sample sheet the numbers
   match under other column names; name series as on GEO (`GM_NCP_sp`, `mCD8T_WT_NCP_sp`) or give `chroms`,
   so chrY is left out where Park et al. left it out.

```bash
condense fetch SRR31626986 SRR31626987 SRR31626988 SRR31626989 -o fastq --threads 4
condense reference hg38_ucsc --bowtie2 2.4.1 --threads 8    # once, ~1.5 h; prints the FASTA and index paths
condense config.yaml --park
```

On a cluster, run the first two in a job (`condense` prints an `sbatch` line to paste) and add `--slurm` to the third.
`config.yaml` for this example (H1 PEG rep 2) is in [details.md](details.md#config-files).

Tables tested identical to GEO's, byte for byte:
`H1_NCP_PEG_2rep` (10 kb), `H1_NCP_sp_2rep` (10 kb), `H1_NCP_sp_1rep_deep` (1 kb and per nucleosome),
`E14_NCP_sp_1rep` and `E14_NCP_sp_2rep` (10 kb).

## Changing one Park setting

Settings you give override `--park`. For example, Park et al.'s processing without multi-mapping reads:

```bash
echo 'reads: {min_mapq: 10}' > settings.yaml
condense samples.tsv --park --settings settings.yaml
```

In a config run, put the setting in the config file. All settings: [details.md](details.md#all-settings).

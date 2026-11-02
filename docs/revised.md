# The revised processing (`--revised`)

By default condense-calc reproduces the processing of Park & Ha, Nature, 2025 (GEO GSE252941; scripts:
[spark159/condense-seq](https://github.com/spark159/condense-seq)): from Park et al.'s reads, the tables are
identical to GEO's ([geo.md](geo.md)). `--revised` (or `park: false` in a settings or config file) changes a
few design choices.

```bash
condense samples.tsv --revised
```

## What `--revised` changes

| what | default (Park et al.) | `--revised` | why |
|---|---|---|---|
| Reads that fit several places in the genome equally well (multi-mapping; `reads: {min_mapq}`) | kept (0: no filter) | dropped: a fragment counts only if its forward-strand read has MAPQ ≥ 10. Some 10 kb bins (mostly centromeres and other repeats) then have no input reads and are NA | Where bowtie2 puts such a read depends on the read names and base qualities, so their counts differ between pipelines. |
| Bins without input reads (`score: {missing_as_nan}`) | dropped | written as NA rows | Every bin table has every bin, so bin tables from different runs line up row for row. |
| Human reference (`-g hg38`) | UCSC hg38.fa (`hg38_ucsc`: 455 sequences, 261 of them alt haplotypes); the index is built on first use (~1.5 h) | NCBI GRCh38 no-alt analysis set, with NCBI's prebuilt index | Alt haplotypes are extra copies of parts of the main chromosomes, so bowtie2 sees reads there as multi-mapping. |
| bowtie2 (`bowtie2`) | 2.4.1 (Park et al.'s Reporting Summary), installed once as a conda env in the refs folder | 2.5.5, the environment's | The environment's version needs no second install. On the H1 libraries tested the two gave the same counts; on E14 rep 1, 2.5.5 differed by one count in 10 bins. |
| FASTQ check (`fastq_check`) | `park`: stop on SRA Lite reads (base qualities all `?`) and on read names not in the form Park et al. aligned (`@<spot> 1:N:0:<index>`) | `lite`: stop on SRA Lite reads only | Read names change where bowtie2 puts reads with several equally good hits; with `--revised` almost all of these have MAPQ < 10 and are dropped. |
| Nucleosome calling: positions with equal counts | in the order Python 2's sort left them | taken in order of position | Park's order is a side effect of Python 2's sort, not a rule. |
| Fragment centre | Park's formula | from the fragment's two ends | Park's formula is one base off for some odd-length fragments and wraps a centre before position 0 into the last bin. bowtie2's end-to-end paired alignment produces neither case, so no output changes. |
| chrY, when `chroms` is not given | dropped for GM (GM12878) and mCD8T, matched by the first part of the series or cell name (`GM_...`, `mCD8T_...`) | kept: all main chromosomes | Park's driver did this by name. `--revised` does not guess from names; give `chroms` to drop chrY. |

Mouse (`-g mm10`) is the ENCODE mm10 no-alt analysis set in both modes: equivalent to the index Park et al. used
(the E14 tables reproduce exactly).

The same in both modes:
- Soluble fractions are used at full precision, exactly as written in the sample sheet or titration file.
- Misspelled or invalid settings stop the run with a message.
- Each run writes `run_info.json`: mode (`park` or `revised`), settings, references, and the bowtie2 version of each library.
- Sliding-window areas of deep series (Park's `Bdata` files) are made only with `bdata: true`. No GEO table uses them.
- Chromosomes of human and mouse genomes are in GEO order (1..N, X, Y), whatever the order in the FASTA.
- The steps and their parameters ([details.md](details.md#steps)), the score formula (with 1 added to every
  count), the c½ fit, no duplicate removal, and at most 10 mismatched or gapped bases per read.

## Changing one setting

Settings you give override both modes. For example, Park et al.'s processing without multi-mapping reads:

```bash
echo 'reads: {min_mapq: 10}' > settings.yaml
condense samples.tsv --settings settings.yaml
```

In a config run, put the setting in the config file. All settings: [details.md](details.md#all-settings).

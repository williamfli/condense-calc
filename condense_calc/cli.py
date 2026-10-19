"""condense-seq reads -> condensability tables.

  condense init FASTQ_DIR      write a sample sheet listing the FASTQ files (fill in conc, fraction)
  condense SHEET               run everything: alignment, scores, c-half fit, QC report
  condense SHEET --slurm       the same, each step as a Slurm job

The other subcommands are single pipeline stages, called by the workflow."""
import argparse
import json
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam
import yaml

from . import counts, fit, fragments, gtab, qc, reference, score, sheet, titration

PKG = Path(__file__).resolve().parent


def label_of(bam):
    """Column label of a BAM, as the original scripts made it (file name up to the first dot)."""
    return Path(bam).name.split('.')[0]


def chrom_sizes(fasta):
    with pysam.FastaFile(str(fasta)) as fa:
        return dict(zip(fa.references, fa.lengths))


def _read(args, chroms):
    return fragments.read_bam(args.bam, chroms, mm_cutoff=args.mm, min_mapq=args.min_mapq, dedup=args.dedup)


def cmd_bincount(args):
    sizes = chrom_sizes(args.fasta)
    chroms = gtab.sort_chroms(args.chroms)
    frags = _read(args, chroms)
    gc = counts.gc_content(args.fasta, chroms, args.win) if args.gc else None
    rows = counts.bincount(frags, sizes, chroms, args.win, gc=gc, tlen=args.tlen)
    header = ['Chromosome', 'Start', 'End', args.label or label_of(args.bam)]
    header += (['GCcontent'] if args.gc else []) + (['Meantlen'] if args.tlen else [])
    gtab.write(args.out, header, rows)


def cmd_coverage(args):
    sizes = chrom_sizes(args.fasta)
    frags = _read(args, [args.chrom])
    f = frags.get(args.chrom)
    cov = counts.coverage(f, sizes[args.chrom]) if f is not None else []
    rows = counts.coverage_rows(args.chrom, cov) if f is not None else []
    gtab.write(args.out, ['Chromosome', 'Position', args.label or label_of(args.bam)], rows)


def cmd_peaks(args):
    frags = _read(args, [args.chrom])
    f = frags.get(args.chrom)
    pk = counts.peaks(f, fast=args.fast) if f is not None else {}
    rows = [[args.chrom, p, pk[p]] for p in sorted(pk)]
    gtab.write(args.out, ['Chromosome', 'Position', args.label or label_of(args.bam)], rows)


def cmd_areas(args):
    """Per-nucleosome (Ncov) and sliding-window (Bdata) coverage areas of one chromosome."""
    sizes = chrom_sizes(args.fasta)
    label = args.label or label_of(args.bam)
    frags = _read(args, [args.chrom])
    f = frags.get(args.chrom)
    size = sizes[args.chrom]
    cov = counts.coverage(f, size) if f is not None else None
    peak_table = gtab.read(args.peaks)
    dyads = [int(r[1]) for r in peak_table.rows if r[0] == args.chrom and float(r[-1]) > 0]
    if args.ncov_out:
        vals = counts.ncpcov(cov, dyads, args.window) if cov is not None else [0] * len(dyads)
        gtab.write(args.ncov_out, ['Chromosome', 'Position', label],
                   [[args.chrom, d, int(v)] for d, v in zip(dyads, vals)])
    if args.bdata_out:
        rows = []
        if cov is not None:
            st, ed, sm = counts.bindata(cov, size, args.bsize, args.bstep)
            rows = [[args.chrom, int(a), int(b), float(v)] for a, b, v in zip(st, ed, sm)]
        gtab.write(args.bdata_out, ['Chromosome', 'Start', 'End', label], rows)


def cmd_score(args):
    if len(args.test) != len(args.control) or len(args.test) != len(args.out):
        sys.exit('--test, --control and --out need the same number of files')
    f = None
    if args.titration:
        f = titration.read(args.titration).frac[args.tnum]
    res = score.score(list(zip(args.test, args.control)), survival_fraction=f,
                      chroms=args.chroms, missing_as_nan=args.missing_as_nan)
    for (header, rows), out in zip(res, args.out):
        gtab.write(out, header, rows)


def cmd_num(args):
    f = titration.read(args.titration).frac[args.tnum]
    res = score.numbers(args.inputs, f, mscale=args.mscale, chroms=args.chroms)
    for (header, rows), out in zip(res, args.out):
        gtab.write(out, header, rows)


def cmd_concat(args):
    header, rows = score.concatenate(args.inputs)
    gtab.write(args.out, header, rows)


def cmd_fit(args):
    tit = titration.read(args.titration)
    header, rows = fit.fit_table(args.input, tit, args.tnums, model=args.model, method=args.method,
                                 min_rsq=args.min_rsq, chroms=args.chroms)
    gtab.write(args.out, header, rows)


def cmd_init(args):
    out = Path(args.out)
    if out.exists() and not args.force:
        sys.exit(f'{out} exists; use --force to overwrite')
    fastqs = sheet.find_fastqs(args.dirs) if args.dirs else []
    if args.dirs and not fastqs:
        sys.exit(f'no read-1 FASTQ files (*_R1*.fastq.gz, *_1.fastq.gz) in {" ".join(args.dirs)}')
    lines = ['# One row per read-1 FASTQ (read 2 is found by name). Rows with the same series and conc are',
             '# lanes of one library. conc: agent concentration, 0 for the input. fraction: soluble fraction',
             '# measured at that concentration (input: 1). Optional columns: rep, agent, fastq2.',
             'series\tconc\tfraction\tfastq']
    if fastqs:
        here = out.resolve().parent
        lines += [f'{args.series}\t\t\t{f.resolve().relative_to(here) if f.resolve().is_relative_to(here) else f.resolve()}'
                  for f in fastqs]
    else:
        lines += [f'{args.series}\t0\t1\treads/input_R1.fastq.gz',
                  f'{args.series}\t0.794\t0.105\treads/sp8_L001_R1.fastq.gz',
                  f'{args.series}\t0.794\t0.105\treads/sp8_L002_R1.fastq.gz']
    out.write_text('\n'.join(lines) + '\n')
    print(f'wrote {out}' + (f' with {len(fastqs)} FASTQ files; fill in conc and fraction' if fastqs else ''))


def cmd_reference(args):
    for g in args.genomes:
        print('\t'.join(reference.ensure(g, args.refs, args.threads)))


def _check(src):
    with tempfile.TemporaryDirectory() as tmp:
        sheet.read_simple(src, 'check', workdir=tmp)


def _workflow_config(args, out):
    """Config for the workflow; a sample sheet becomes condense.yaml in the output folder."""
    src = Path(args.sheet).resolve()
    if src.suffix in ('.yaml', '.yml'):
        return src
    _check(src)                                                      # fail early, before any job starts
    name, fasta, index = reference.ensure(args.genome, args.refs, args.cores)
    cfg = yaml.safe_load(Path(args.settings).read_text()) if args.settings else {}
    cfg.update(sheet=str(src), genome=name, genomes={name: {'fasta': fasta, 'bowtie2_index': index}},
               nucleosomes=args.nucleosomes)
    path = out / 'condense.yaml'
    text = yaml.safe_dump(cfg, sort_keys=False)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)
    return path


def _submit_controller(args, out):
    """Run this command as a small Slurm job that submits and watches the steps."""
    argv = [sys.executable, '-m', 'condense_calc', *sys.argv[1:]]
    # controller: 1 core, < 1 GB; it runs as long as the whole workflow
    cmd = ['sbatch', '--parsable', '-J', 'condense', '-c', '1', '--mem', '2G', '-t', args.time,
           '-o', str(out / 'condense_%j.log')]
    cmd += ['-A', args.account] if args.account else []
    cmd += ['-p', args.partition] if args.partition else []
    cmd += ['--wrap', shlex.join(argv)]
    job = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()
    print(f'Submitted job {job}. Progress: {out}/condense_{job}.log\n'
          f'Results will be in {out}/tables and {out}/qc/report.html')


def cmd_run(args):
    src = Path(args.sheet).resolve()
    if not src.exists():
        sys.exit(f'{args.sheet} not found')
    out = Path(args.outdir).resolve() if args.outdir else src.parent / 'results'
    out.mkdir(parents=True, exist_ok=True)
    if args.slurm and not args.dry_run and 'SLURM_JOB_ID' not in os.environ:
        if src.suffix not in ('.yaml', '.yml'):
            try:
                _check(src)
            except sheet.SheetError as e:
                sys.exit(str(e))
        return _submit_controller(args, out)
    try:
        cfg = _workflow_config(args, out)
    except sheet.SheetError as e:
        sys.exit(str(e))
    cmd = ['snakemake', '-s', str(PKG / 'workflow' / 'Snakefile'), '--configfile', str(cfg),
           '--directory', str(out), '--rerun-incomplete']
    if args.qc:
        cmd += ['--config', f'qc={args.qc}']
    if args.slurm:
        res = ['mem_mb=4000', 'runtime=120']
        res += [f'slurm_account={args.account}'] if args.account else []
        res += [f'slurm_partition={args.partition}'] if args.partition else []
        # --cores is the thread budget, not the controller's cores: without it Snakemake takes the CPUs of
        # the (1-core) controller job and caps every rule's threads at 1, so bowtie2 jobs get 1 CPU each.
        # --jobs limits how many jobs run at once; this only has to be large enough not to bind.
        cmd += ['--executor', 'slurm', '--jobs', str(args.jobs), '--cores', str(args.jobs * 64),
                '--retries', '1', '--latency-wait', '60', '--default-resources', *res]
    else:
        cmd += ['--cores', str(args.cores)]
    if args.dry_run:
        cmd.append('--dry-run')
    if args.quiet:
        cmd += ['--quiet', 'rules']
    cmd += args.snakemake_args
    print('+ ' + ' '.join(cmd), flush=True)
    rc = subprocess.call(cmd)
    flags = out / 'qc' / 'flags.tsv'
    if rc == 0 and not args.dry_run:
        msg = f'\nTables: {out / "tables"}'
        if flags.exists():
            n = sum(1 for _ in open(flags)) - 1
            msg += f'\nQC report: {out / "qc" / "report.html"} ({n} flag{"s" if n != 1 else ""}; see {flags})'
        print(msg)
    sys.exit(rc)


def cmd_qc_sample(args):
    d = qc.sample_stats(args.bam, args.bowtie2_log, mm=args.mm, min_mapq=args.min_mapq,
                        threads=args.threads, verbose=args.verbose)
    qc._write_json(args.out, d)


def cmd_qc_series(args):
    d = qc.series_stats(args.series, args.samples, args.bins, args.table, args.titration, fit=args.fit,
                        peaks=args.peaks or (), fasta=args.fasta, verbose=args.verbose)
    qc._write_json(args.out, d)
    qc.series_figure(d, [json.loads(Path(p).read_text()) for p in args.samples], args.png)


def cmd_qc_replicates(args):
    d = qc.replicate_stats(args.group, {t[0]: t[1:] for t in args.tables})
    qc._write_json(args.out, d)
    qc.replicate_figure(d, args.png)


def cmd_qc_report(args):
    fl = qc.report(args.samples, args.series, args.replicates or [], args.thresholds or [], args.outdir)
    for f in fl:
        v = qc._fmt(f['value']) if not isinstance(f['value'], str) else f['value']
        print(f"QC {f['level']}: {f['name']}: {f['metric']} = {v} ({f['threshold']})")


def main(argv=None):
    p = argparse.ArgumentParser(prog='condense', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    s = sub.add_parser('init', help='write a sample sheet listing the FASTQ files in some folders')
    s.add_argument('dirs', nargs='*', help='folders with FASTQ files (none: write an example sheet)')
    s.add_argument('-o', '--out', default='samples.tsv', help='sheet to write (default samples.tsv)')
    s.add_argument('--series', default='series1', help='series name to fill in (default series1)')
    s.add_argument('--force', action='store_true', help='overwrite an existing sheet')
    s.set_defaults(func=cmd_init)

    s = sub.add_parser('reference', help='download a genome (FASTA + bowtie2 index); run does this when needed')
    s.add_argument('genomes', nargs='+', help='hg38, mm10, or a FASTA file')
    s.add_argument('--refs', help='where genomes are kept (default $CONDENSE_REFS, else ~/.cache/condense-calc)')
    s.add_argument('--threads', type=int, default=1, help='for building an index (mm10, own FASTA)')
    s.set_defaults(func=cmd_reference)

    s = sub.add_parser('run', help='run everything (the default command)',
                       description='Run everything. Arguments after -- go to snakemake.')
    s.add_argument('sheet', help='sample sheet (condense init), or a config.yaml for multi-genome runs')
    s.add_argument('-o', '--outdir', help='output folder (default: results/ next to the sheet)')
    s.add_argument('-g', '--genome', default='hg38', help='hg38 (default), mm10, or a FASTA file')
    s.add_argument('--refs', help='where genomes are kept (default $CONDENSE_REFS, else ~/.cache/condense-calc)')
    s.add_argument('--nucleosomes', action='store_true',
                   help='also score every nucleosome (for deeply sequenced series; slower)')
    s.add_argument('--qc', choices=['off', 'default', 'verbose'], help='QC level (default: default)')
    s.add_argument('--settings', help='YAML file with further settings (see README)')
    s.add_argument('--cores', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK') or os.cpu_count() or 4),
                   help='local run: cores to use (default: all)')
    s.add_argument('--slurm', action='store_true', help='submit each step as a Slurm job')
    s.add_argument('--account', help='Slurm account')
    s.add_argument('--partition', help='Slurm partition(s), e.g. shared or shared,sapphire')
    s.add_argument('--time', default='3-00:00', help='Slurm: time limit of the controlling job (default 3 days)')
    s.add_argument('--jobs', type=int, default=100, help='Slurm: max jobs at once (default 100)')
    s.add_argument('-q', '--quiet', action='store_true', help='less progress output')
    s.add_argument('-n', '--dry-run', action='store_true', help='show what would run')
    s.set_defaults(func=cmd_run)

    def reads(sp):
        sp.add_argument('--bam', required=True)
        sp.add_argument('--label', help='column name (default: BAM file name up to the first dot)')
        sp.add_argument('--mm', type=int, default=10, help='max edit distance NM (default 10, as GEO)')
        sp.add_argument('--min-mapq', type=int, default=0,
                        help='MAPQ filter; 0 = none, as the original (departure if > 0)')
        sp.add_argument('--dedup', action='store_true',
                        help='drop duplicate fragments (departure; the original keeps them)')

    s = sub.add_parser('bincount', help='reads per genomic bin (bincount_ver3.py)')
    reads(s)
    s.add_argument('--fasta', required=True)
    s.add_argument('--chroms', nargs='+', required=True)
    s.add_argument('--win', type=int, required=True)
    s.add_argument('--no-gc', dest='gc', action='store_false')
    s.add_argument('--no-tlen', dest='tlen', action='store_false')
    s.add_argument('-o', '--out', required=True)
    s.set_defaults(func=cmd_bincount)

    s = sub.add_parser('coverage', help='fragment coverage of one chromosome (coverage_ver3.py --skip)')
    reads(s)
    s.add_argument('--fasta', required=True)
    s.add_argument('--chrom', required=True)
    s.add_argument('-o', '--out', required=True)
    s.set_defaults(func=cmd_coverage)

    s = sub.add_parser('peaks', help='nucleosome dyads of one chromosome (NCPpeak_ver3.py)')
    reads(s)
    s.add_argument('--chrom', required=True)
    s.add_argument('--fast', action='store_true', help='break score ties by position (not the original order)')
    s.add_argument('-o', '--out', required=True)
    s.set_defaults(func=cmd_peaks)

    s = sub.add_parser('areas', help='coverage area per dyad (NCPcov_ver3.py) and per sliding window (bindata_ver3.py)')
    reads(s)
    s.add_argument('--fasta', required=True)
    s.add_argument('--chrom', required=True)
    s.add_argument('--peaks', required=True, help='peak gtab of the input sample')
    s.add_argument('--window', type=int, default=counts.NCP_WINDOW)
    s.add_argument('--bsize', type=int, default=counts.SLIDING[0])
    s.add_argument('--bstep', type=int, default=counts.SLIDING[1])
    s.add_argument('--ncov-out')
    s.add_argument('--bdata-out')
    s.set_defaults(func=cmd_areas)

    s = sub.add_parser('score', help='condensability scores of one titration point (NCPscore_ver4.py --numc)')
    s.add_argument('--test', nargs='+', required=True)
    s.add_argument('--control', nargs='+', required=True)
    s.add_argument('--titration', help='titration table; omit for no survival-fraction correction')
    s.add_argument('--tnum', type=int)
    s.add_argument('--chroms', nargs='+')
    s.add_argument('--missing-as-nan', action='store_true',
                   help='write bins without input reads as NA instead of dropping them (departure)')
    s.add_argument('-o', '--out', nargs='+', required=True)
    s.set_defaults(func=cmd_score)

    s = sub.add_parser('num', help='nucleosome numbers of one titration point (NCPnum_ver2.py)')
    s.add_argument('--inputs', nargs='+', required=True)
    s.add_argument('--titration', required=True)
    s.add_argument('--tnum', type=int, required=True)
    s.add_argument('--mscale', type=float, default=1)
    s.add_argument('--chroms', nargs='+')
    s.add_argument('-o', '--out', nargs='+', required=True)
    s.set_defaults(func=cmd_num)

    s = sub.add_parser('concat', help='outer join of gtabs (concat_gtab.py)')
    s.add_argument('--inputs', nargs='+', required=True)
    s.add_argument('-o', '--out', required=True)
    s.set_defaults(func=cmd_concat)

    s = sub.add_parser('fit', help='c-half logistic fit per row (logistic_fit_ver3.py)')
    s.add_argument('--input', required=True)
    s.add_argument('--titration', required=True)
    s.add_argument('--tnums', type=int, nargs='+', required=True)
    s.add_argument('--model', default='sigmoid', choices=sorted(fit.MODELS))
    s.add_argument('--method', default='evolution', choices=['evolution', 'curve_fit'])
    s.add_argument('--min-rsq', type=float, default=0.5)
    s.add_argument('--chroms', nargs='+')
    s.add_argument('-o', '--out', required=True)
    s.set_defaults(func=cmd_fit)

    s = sub.add_parser('qc-sample', help='QC: alignment, read filters and fragment lengths of one BAM')
    s.add_argument('--bam', required=True)
    s.add_argument('--bowtie2-log', required=True)
    s.add_argument('--mm', type=int, default=10)
    s.add_argument('--min-mapq', type=int, default=0)
    s.add_argument('--threads', type=int, default=1)
    s.add_argument('--verbose', action='store_true')
    s.add_argument('-o', '--out', required=True)
    s.set_defaults(func=cmd_qc_sample)

    s = sub.add_parser('qc-series', help='QC of one titration series')
    s.add_argument('--series', required=True)
    s.add_argument('--samples', nargs='+', required=True, help='qc-sample JSON files')
    s.add_argument('--bins', nargs='+', required=True, help='bin-count tables of all points, input included')
    s.add_argument('--table', required=True, help='binned score table')
    s.add_argument('--titration', required=True)
    s.add_argument('--fit')
    s.add_argument('--peaks', nargs='*')
    s.add_argument('--fasta')
    s.add_argument('--verbose', action='store_true')
    s.add_argument('-o', '--out', required=True)
    s.add_argument('--png', required=True)
    s.set_defaults(func=cmd_qc_series)

    s = sub.add_parser('qc-replicates', help='QC: agreement between replicate series')
    s.add_argument('--group', required=True)
    s.add_argument('--tables', nargs='+', action='append', required=True, metavar='KB_THEN_FILES',
                   help='a bin-size label then one score table per replicate; repeat per bin size')
    s.add_argument('-o', '--out', required=True)
    s.add_argument('--png', required=True)
    s.set_defaults(func=cmd_qc_replicates)

    s = sub.add_parser('qc-report', help='QC: summary, flags and HTML report')
    s.add_argument('--samples', nargs='+', required=True)
    s.add_argument('--series', nargs='+', required=True)
    s.add_argument('--replicates', nargs='*')
    s.add_argument('--thresholds', nargs='*', help='overrides, e.g. min_replicate_rho=0.7')
    s.add_argument('--outdir', required=True)
    s.set_defaults(func=cmd_qc_report)

    argv = list(sys.argv[1:] if argv is None else argv)
    extra = argv[argv.index('--') + 1:] if '--' in argv else []     # `condense-calc run ... -- <snakemake args>`
    argv = argv[:argv.index('--')] if '--' in argv else argv
    if argv and not argv[0].startswith('-') and argv[0] not in sub.choices:
        argv = ['run'] + argv                                         # `condense SHEET` = `condense run SHEET`
    args = p.parse_args(argv)
    args.snakemake_args = extra
    args.func(args)


if __name__ == '__main__':
    main()

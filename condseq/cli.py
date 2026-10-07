"""Command line: `python -m condseq <stage> ...`. One subcommand per pipeline stage; the
Snakemake workflow (workflow/Snakefile) calls these. Every stage defaults to the settings of
the Park 2025 GEO run; the opt-in departures are --min-mapq, --dedup and --missing-as-nan."""
import argparse
import sys
from pathlib import Path

import pysam

from . import counts, fit, fragments, gtab, score, titration


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


def main(argv=None):
    p = argparse.ArgumentParser(prog='condseq', description=__doc__)
    sub = p.add_subparsers(dest='cmd', required=True)

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

    args = p.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()

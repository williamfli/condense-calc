"""condense-seq reads -> condensability tables.

  condense init FASTQ_DIR      write a sample sheet listing the FASTQ files (fill in conc, fraction)
  condense SHEET               run everything: alignment, scores, c-half fit, QC report
  condense SHEET --slurm       the same, each step as a Slurm job

The other subcommands are single pipeline stages, called by the workflow."""
import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam
import yaml

from . import counts, fastq, fit, fragments, gtab, qc, reference, score, settings, sheet, titration

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
    rows = counts.bincount(frags, sizes, chroms, args.win, gc=gc, tlen=args.tlen, park=args.park)
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
    pk = counts.peaks(f, park=args.park) if f is not None else {}
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
    if not args.genomes and not args.bowtie2:
        sys.exit('nothing to do: name a genome, or give --bowtie2 VERSION')
    for g in args.genomes:
        print('\t'.join(reference.ensure(g, args.refs, args.threads)))
    if args.bowtie2 and reference.missing_bowtie2(args.bowtie2, args.refs) and _on_login_node():
        sys.exit(f'Installing bowtie2 {args.bowtie2} (conda) is not light work for a login node; run it in a job, '
                 f"e.g. sbatch -c 1 --mem 8G -t 1:00:00 --wrap 'condense reference --bowtie2 {args.bowtie2}', "
                 'or let condense ... --slurm do it.')
    if args.bowtie2:
        d = reference.ensure_bowtie2(args.bowtie2, args.refs)
        print(f'bowtie2 {reference.bowtie2_version(d)}: {d}' if d else 'bowtie2: the one on PATH')


def cmd_check_fastq(args):
    res = fastq.check(args.r1.split(','), args.r2.split(','), args.mode)
    if args.out:
        Path(args.out).write_text(json.dumps({'mode': args.mode, 'files': res}, indent=1) + '\n')
    if fastq.problems(res):
        sys.exit('\n'.join(fastq.problems(res)))


def _check(src, fastq_check='off', out=None):
    """Check a sample sheet (or a config's samples), and the first reads of the FASTQ files of libraries not
    yet aligned in out; raises SheetError. The align step checks the FASTQs again."""
    if src.suffix in ('.yaml', '.yml'):
        cfg = yaml.safe_load(src.read_text()) or {}
        if 'samples' not in cfg:
            return
        libs = sheet.read_config_sheet(cfg['samples'], cfg.get('titration_dir', 'titration'),
                                       lambda p: str((src.parent / Path(p).expanduser()).resolve()))
    else:
        with tempfile.TemporaryDirectory() as tmp:
            libs = sheet.read_simple(src, 'check', workdir=tmp, note=lambda m: print(m, file=sys.stderr))
    todo = [r for name, r in libs.iterrows() if not (out and (Path(out) / 'bam' / f'{name}.bam').exists())]
    found = lambda fs: [f for f in fs.split(',') if Path(f).exists()]     # a missing one: Snakemake says so
    probs = [p for r in todo for p in fastq.problems(fastq.check(found(r['fastq1']), found(r['fastq2']), fastq_check))]
    if probs:
        raise sheet.SheetError(f'FASTQ check (fastq_check: {fastq_check}):\n  ' + '\n  '.join(probs))


def _user_settings(args, src):
    """The user's settings: a config file's own keys, or the --settings file of a sample-sheet run; plus --park."""
    path = src if src.suffix in ('.yaml', '.yml') else (Path(args.settings) if args.settings else None)
    cfg = (yaml.safe_load(path.read_text()) or {}) if path else {}
    if args.park:
        cfg['park'] = True
    return cfg


def _todo(args, S, is_sheet):
    """What has to be prepared before the workflow can start: the genome of a sheet run, a bowtie2 version."""
    todo = [reference.missing(args.genome, args.refs)] if is_sheet else []
    return [t for t in todo + [reference.missing_bowtie2(S['bowtie2'], args.refs)] if t]


def _workflow_config(args, out):
    """Config for the workflow; a sample sheet becomes condense.yaml in the output folder."""
    src = Path(args.sheet).resolve()
    if src.suffix in ('.yaml', '.yml'):
        return src
    name, fasta, index = reference.ensure(args.genome, args.refs, _local_cores(args))
    cfg = _user_settings(args, src)
    cfg.update(sheet=str(src), genome=name, genomes={name: {'fasta': fasta, 'bowtie2_index': index}},
               refs=str(reference.refs_dir(args.refs)))
    if args.nucleosomes:                                             # else keep the --settings file's value
        cfg['nucleosomes'] = True
    path = out / 'condense.yaml'
    text = yaml.safe_dump(cfg, sort_keys=False)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)
    return path


CONTROLLER = 'CONDENSE_CONTROLLER'     # set in the controller job, which runs the workflow itself


def _local_cores(args):
    """--cores, else the CPUs this process may use (in a Slurm job: its allocation, not the node)."""
    return args.cores or len(os.sched_getaffinity(0))


def _local_mem_mb(cores):
    """Memory (MB) for a local run: the Slurm allocation if inside one, else 90% of the machine."""
    env = os.environ
    if env.get('SLURM_MEM_PER_NODE'):
        return int(env['SLURM_MEM_PER_NODE'])
    if env.get('SLURM_MEM_PER_CPU'):
        return int(env['SLURM_MEM_PER_CPU']) * cores
    try:
        return int(os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') * 0.9 / 1e6)
    except (ValueError, OSError, AttributeError):
        return None


def _on_login_node():
    """A cluster (sbatch exists) but not inside a job: e.g. a Cannon login node, meant for light work."""
    return 'SLURM_JOB_ID' not in os.environ and shutil.which('sbatch') is not None


def _sbatch(args, cmd):
    cmd = cmd[:1] + ['--parsable'] + (['-A', args.account] if args.account else []) + \
        (['-p', args.partition] if args.partition else []) + \
        (['--mail-type=FAIL,TIME_LIMIT_80', f'--mail-user={args.mail}'] if args.mail else []) + cmd[1:]
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()


def _submit_controller(args, out, S):
    """Run this command as a small Slurm job that submits and watches the steps. A genome that still has
    to be downloaded or indexed, or a bowtie2 version to install, is prepared first in its own job (8 cores),
    which the controller waits for."""
    dep = []
    is_sheet = Path(args.sheet).suffix not in ('.yaml', '.yml')
    todo = _todo(args, S, is_sheet)
    if todo:
        ref = [sys.executable, '-m', 'condense_calc', 'reference', '--threads', '8', '--bowtie2', S['bowtie2']]
        if is_sheet:
            ref.append(str(reference.locate(args.genome)[1]) if args.genome not in reference.PRESETS else args.genome)
        ref += ['--refs', args.refs] if args.refs else []
        # bowtie2-build of a mammalian genome: ~1-2 h on 8 cores, < 16 GB
        rjob = _sbatch(args, ['sbatch', '-J', 'condense-reference', '-c', '8', '--mem', '32G', '-t', '0-12:00',
                              '-o', str(out / 'condense_reference_%j.log'), '--wrap', shlex.join(ref)])
        print(f'Submitted job {rjob} to {" and ".join(todo)}; the run starts when it is done.')
        dep = [f'--dependency=afterok:{rjob}', '--kill-on-invalid-dep=yes']
    argv = [sys.executable, '-m', 'condense_calc', *sys.argv[1:]]
    # controller: 1 core, < 1 GB; it runs as long as the whole workflow
    job = _sbatch(args, ['sbatch', '-J', 'condense', '-c', '1', '--mem', '2G', '-t', args.time, *dep,
                         f'--export=ALL,{CONTROLLER}=1', '-o', str(out / 'condense_%j.log'),
                         '--wrap', shlex.join(argv)])
    print(f'Submitted job {job}. Progress: {out}/condense_{job}.log\n'
          f'Results will be in {out}/tables and {out}/qc/report.html')


def cmd_run(args):
    src = Path(args.sheet).resolve()
    if not src.exists():
        sys.exit(f'{args.sheet} not found')
    is_sheet = src.suffix not in ('.yaml', '.yml')
    try:
        S = settings.resolve(_user_settings(args, src))             # unknown keys: stop before anything starts
    except (settings.SettingsError, OSError, yaml.YAMLError) as e:
        sys.exit(f'settings: {e}')
    if S['park'] and is_sheet and args.genome in settings.PARK_GENOMES:
        args.genome = settings.PARK_GENOMES[args.genome]
        print(f'Park mode: genome {args.genome} ({reference.PRESETS[args.genome]["name"]})', flush=True)
    out = Path(args.outdir).resolve() if args.outdir else src.parent / 'results'
    out.mkdir(parents=True, exist_ok=True)
    try:
        _check(src, S['fastq_check'], out)                          # fail early, before any job starts
    except sheet.SheetError as e:
        sys.exit(str(e))
    # --slurm submits a controller job, also from an interactive (salloc / Open OnDemand) session: run in
    # the foreground there, the workflow would die with the session
    if args.slurm and not args.dry_run and os.environ.get(CONTROLLER) != '1':
        return _submit_controller(args, out, S)
    if not args.slurm and not args.dry_run and args.cores is None and _on_login_node():
        sys.exit('This looks like a cluster login node (sbatch is available, but this is not a Slurm job).\n'
                 'Login nodes are for light work: add --slurm to run each step as a Slurm job, or run inside\n'
                 'a job (salloc / sbatch). To run here anyway, give --cores N.')
    todo = _todo(args, S, is_sheet)
    if args.dry_run and todo:
        print(f'{"Sample sheet OK. " if is_sheet else ""}Dry run stops here: the run would first '
              f'{" and ".join(todo)}.\nRun without -n to do that (with --slurm it runs as its own job).')
        sys.exit(0)
    # a bowtie2 version is installed by the preparation job (--slurm), or here in a local run inside a job
    bt2 = reference.missing_bowtie2(S['bowtie2'], args.refs)
    if bt2 and (os.environ.get(CONTROLLER) == '1' or _on_login_node()):
        sys.exit(f'First {bt2}: condense reference --bowtie2 {S["bowtie2"]}, in a Slurm job')
    reference.ensure_bowtie2(S['bowtie2'], args.refs)
    try:
        cfg = _workflow_config(args, out)
    except sheet.SheetError as e:
        sys.exit(str(e))
    cmd = ['snakemake', '-s', str(PKG / 'workflow' / 'Snakefile'), '--configfile', str(cfg),
           '--directory', str(out), '--rerun-incomplete']
    conf = [f'qc={args.qc}'] if args.qc else []
    conf += ['park=true'] if args.park and not is_sheet else []     # a sheet run has it in condense.yaml
    conf += [f'refs={reference.refs_dir(args.refs)}'] if args.refs and not is_sheet else []
    cmd += ['--config', *conf] if conf else []
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
        cores = _local_cores(args)
        cmd += ['--cores', str(cores)]
        mem = _local_mem_mb(cores)
        # steps request memory (mem_mb); without a total, Snakemake would start as many as there are
        # cores and could run the machine (or the job's allocation) out of memory
        cmd += ['--resources', f'mem_mb={mem}'] if mem else []
    if args.dry_run:
        cmd.append('--dry-run')
    if args.quiet:
        cmd += ['--quiet', 'rules']
    cmd += args.snakemake_args
    print('+ ' + ' '.join(cmd), flush=True)
    # the tools (snakemake, bowtie2, samtools) of this Python's environment, also when it is not activated
    env = dict(os.environ, PATH=f'{Path(sys.executable).parent}{os.pathsep}{os.environ.get("PATH", "")}')
    rc = subprocess.call(cmd, env=env)
    flags = out / 'qc' / 'flags.tsv'
    if rc == 0 and not args.dry_run:
        msg = f'\nTables: {out / "tables"}\nRun record: {out / "run_info.json"}'
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


# stage option -> (settings section, key)
STAGE_DEFAULTS = {'mm': ('reads', 'mm'), 'min_mapq': ('reads', 'min_mapq'), 'dedup': ('reads', 'dedup'),
                  'missing_as_nan': ('score', 'missing_as_nan'), 'mode': ('fastq_check', None),
                  'model': ('fit', 'model'), 'method': ('fit', 'method'), 'min_rsq': ('fit', 'min_rsq')}


def _stage_defaults(args):
    """Stage options not given on the command line take the settings' value (with --park, the Park preset's)."""
    S = settings.resolve(park=getattr(args, 'park', False))
    for dest, (sec, key) in STAGE_DEFAULTS.items():
        if getattr(args, dest, 0) is None:
            setattr(args, dest, S[sec] if key is None else S[sec][key])


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
    s.add_argument('genomes', nargs='*', help=f'{", ".join(reference.PRESETS)}, or a FASTA file')
    s.add_argument('--refs', help='where genomes are kept (default $CONDENSE_REFS, else ~/.cache/condense-calc)')
    s.add_argument('--threads', type=int, default=1, help='for building an index (hg38_ucsc, own FASTA)')
    s.add_argument('--bowtie2', metavar='VERSION',
                   help='also install this bowtie2 version (a conda env in REFS/envs), e.g. 2.4.1 for --park')
    s.set_defaults(func=cmd_reference)

    s = sub.add_parser('run', help='run everything (the default command)',
                       description='Run everything. Arguments after -- go to snakemake.')
    s.add_argument('sheet', help='sample sheet (condense init), or a config.yaml for multi-genome runs')
    s.add_argument('-o', '--outdir', help='output folder (default: results/ next to the sheet)')
    s.add_argument('-g', '--genome', default='hg38',
                   help='hg38 (default), mm10, hg38_ucsc, or a FASTA file; with --park, hg38 means hg38_ucsc')
    s.add_argument('--refs', help='where genomes are kept (default $CONDENSE_REFS, else ~/.cache/condense-calc)')
    s.add_argument('--nucleosomes', action='store_true',
                   help='also score every nucleosome (for deeply sequenced series; slower)')
    s.add_argument('--qc', choices=['off', 'default', 'verbose'], help='QC level (default: default)')
    s.add_argument('--settings', help='YAML file with further settings (see README)')
    s.add_argument('--park', action='store_true',
                   help='reproduce the GEO tables of Park et al. 2025: their processing (settings `park: true`)')
    s.add_argument('--cores', type=int,
                   help='local run: cores to use (default: the CPUs this process may use, e.g. a job\'s allocation)')
    s.add_argument('--slurm', action='store_true', help='submit each step as a Slurm job')
    s.add_argument('--account', help='Slurm account')
    s.add_argument('--partition', help='Slurm partition(s), e.g. shared or shared,sapphire')
    s.add_argument('--time', default='3-00:00', help='Slurm: time limit of the controlling job (default 3 days)')
    s.add_argument('--jobs', type=int, default=100, help='Slurm: max jobs at once (default 100)')
    s.add_argument('--mail', metavar='ADDRESS',
                   help='Slurm: email ADDRESS if the run fails or reaches 80%% of its time limit')
    s.add_argument('-q', '--quiet', action='store_true', help='less progress output')
    s.add_argument('-n', '--dry-run', action='store_true', help='show what would run')
    s.set_defaults(func=cmd_run)

    # stage options left out take their value from settings.py (with --park: the Park preset), see _stage_defaults
    def park(sp):
        sp.add_argument('--park', action='store_true',
                        help="Park et al. 2025's processing: the Park preset's defaults and Park's behaviour")

    def reads(sp):
        sp.add_argument('--bam', required=True)
        sp.add_argument('--label', help='column name (default: BAM file name up to the first dot)')
        sp.add_argument('--mm', type=int, help=f"max edit distance NM ({settings.describe('reads', 'mm')})")
        sp.add_argument('--min-mapq', type=int,
                        help=f"minimum MAPQ, 0 = none ({settings.describe('reads', 'min_mapq')})")
        sp.add_argument('--dedup', action=argparse.BooleanOptionalAction,
                        help=f"drop duplicate fragments ({settings.describe('reads', 'dedup')})")
        park(sp)

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
    s.add_argument('--missing-as-nan', action=argparse.BooleanOptionalAction,
                   help='write bins without input reads as NA rows instead of dropping them '
                        f"({settings.describe('score', 'missing_as_nan')})")
    park(s)
    s.add_argument('-o', '--out', nargs='+', required=True)
    s.set_defaults(func=cmd_score)

    s = sub.add_parser('num', help='nucleosome numbers of one titration point (NCPnum_ver2.py)')
    s.add_argument('--inputs', nargs='+', required=True)
    s.add_argument('--titration', required=True)
    s.add_argument('--tnum', type=int, required=True)
    s.add_argument('--mscale', type=float, default=1)
    s.add_argument('--chroms', nargs='+')
    park(s)
    s.add_argument('-o', '--out', nargs='+', required=True)
    s.set_defaults(func=cmd_num)

    s = sub.add_parser('concat', help='outer join of gtabs (concat_gtab.py)')
    s.add_argument('--inputs', nargs='+', required=True)
    park(s)
    s.add_argument('-o', '--out', required=True)
    s.set_defaults(func=cmd_concat)

    s = sub.add_parser('fit', help='c-half logistic fit per row (logistic_fit_ver3.py)')
    s.add_argument('--input', required=True)
    s.add_argument('--titration', required=True)
    s.add_argument('--tnums', type=int, nargs='+', required=True)
    s.add_argument('--model', choices=sorted(fit.MODELS), help=settings.describe('fit', 'model'))
    s.add_argument('--method', choices=['evolution', 'curve_fit'], help=settings.describe('fit', 'method'))
    s.add_argument('--min-rsq', type=float, help=settings.describe('fit', 'min_rsq'))
    s.add_argument('--chroms', nargs='+')
    park(s)
    s.add_argument('-o', '--out', required=True)
    s.set_defaults(func=cmd_fit)

    s = sub.add_parser('check-fastq', help='check the first reads of FASTQ files before alignment (fastq.py)')
    s.add_argument('-1', dest='r1', required=True, help='read-1 files, comma-separated')
    s.add_argument('-2', dest='r2', required=True, help='read-2 files, comma-separated')
    s.add_argument('--mode', choices=fastq.MODES, help=settings.describe('fastq_check'))
    park(s)
    s.add_argument('-o', '--out', help='JSON with what was seen in each file')
    s.set_defaults(func=cmd_check_fastq)

    s = sub.add_parser('qc-sample', help='QC: alignment, read filters and fragment lengths of one BAM')
    s.add_argument('--bam', required=True)
    s.add_argument('--bowtie2-log', required=True)
    s.add_argument('--mm', type=int, help=settings.describe('reads', 'mm'))
    s.add_argument('--min-mapq', type=int, help=settings.describe('reads', 'min_mapq'))
    park(s)
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
    _stage_defaults(args)
    args.func(args)


if __name__ == '__main__':
    main()

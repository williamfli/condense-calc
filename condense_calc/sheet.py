"""Sample sheets -> one table of libraries for the workflow.

Two formats are read:

* the sample sheet (`condense init`): one row per FASTQ file, with the series, the agent
  concentration and the soluble fraction. Rows with the same series, rep and conc are lanes of one
  library. The titration table of each series is written from the sheet.
* the GEO-style config: `samples` (cell, sample, agent, tnum, rep, depth, fastq1, fastq2) plus one
  titration file per cell/sample/agent. File and column names then follow the GEO tables.

Both give a table indexed by library name with columns series, group (replicate group), cell,
tnum, deep, agent, fastq1, fastq2 (comma-separated lanes) and titration (path).
"""
import glob
import re
from pathlib import Path

import pandas as pd

NAME_RE = re.compile(r'[A-Za-z0-9_-]+')
R1_PATTERNS = [(r'_R1(?=[._])', '_R2'), (r'_1(?=\.f(ast)?q)', '_2'), (r'\.R1(?=\.)', '.R2'), (r'_r1(?=[._])', '_r2')]


class SheetError(ValueError):
    pass


def mate_of(r1):
    """Read-2 file name from a read-1 file name (…_R1_001.fastq.gz, …_1.fastq.gz, …)."""
    name = Path(r1).name
    for pat, rep in R1_PATTERNS:
        new, n = re.subn(pat, rep, name, count=1)
        if n:
            return str(Path(r1).with_name(new))
    raise SheetError(f'{r1}: cannot tell which file is read 2; add a fastq2 column')


def _expand(pattern, base):
    p = Path(pattern).expanduser()
    p = p if p.is_absolute() else base / p
    if glob.has_magic(str(p)):
        hits = sorted(glob.glob(str(p)))
        if not hits:
            raise SheetError(f'no file matches {pattern}')
        return [str(Path(h).resolve()) for h in hits]
    return [str(p.resolve())]


def _sig_digits(text):
    """Significant digits of a number as written ('0.0977444' -> 6, '0.400' -> 3)."""
    return len(text.strip().lower().lstrip('+-').split('e')[0].replace('.', '').lstrip('0'))


def read_simple(path, cell, nucleosomes=False, workdir='.', note=None):
    """The sample sheet. Writes workdir/titration/<series>.tsv (only when its content changes), with conc and
    fraction exactly as written in the sheet. note: called with a message if fractions look rounded."""
    path = Path(path).resolve()
    base = path.parent
    df = pd.read_csv(path, sep='\t', dtype=str, comment='#', skip_blank_lines=True).fillna('')
    df.columns = [c.strip().lower() for c in df.columns]
    for col in ('fastq', 'conc'):
        if col not in df.columns:
            raise SheetError(f'{path}: needs a column {col!r} (columns: series, conc, fraction, fastq)')
    if 'series' not in df.columns:
        df['series'] = path.stem
    has_rep = 'rep' in df.columns and (df['rep'] != '').any()
    for col, default in [('rep', '1'), ('fraction', ''), ('fastq2', ''), ('agent', ''), ('point', '')]:
        if col not in df.columns:
            df[col] = default
    df = df.apply(lambda c: c.str.strip())
    df = df[df['fastq'] != '']

    problems = []
    rows = []
    for _, r in df.iterrows():
        line = r['fastq']
        if not NAME_RE.fullmatch(r['series']):
            problems.append(f'{line}: series {r["series"]!r} may only contain letters, digits, _ and -')
        try:
            conc = float(r['conc'])
        except ValueError:
            problems.append(f'{line}: conc {r["conc"]!r} is not a number (0 for the input)')
            continue
        frac_text = r['fraction']
        if conc == 0:
            frac_text = frac_text or '1'
        try:
            frac = float(frac_text)
            if not 0 < frac <= 1:
                raise ValueError
        except ValueError:
            problems.append(f'{line}: fraction {r["fraction"]!r} must be in (0, 1]; it is the soluble '
                            'fraction measured at this concentration')
            continue
        r1s = _expand(r['fastq'], base)
        r2s = _expand(r['fastq2'], base) if r['fastq2'] else [mate_of(f) for f in r1s]
        if len(r1s) != len(r2s):
            problems.append(f'{line}: {len(r1s)} read-1 files but {len(r2s)} read-2 files')
        for f in r1s + r2s:
            if not Path(f).exists():
                problems.append(f'{line}: {f} not found')
        rows.append(dict(series=r['series'], rep=r['rep'] or '1', conc=conc, fraction=frac,
                         conc_text=r['conc'], frac_text=frac_text,
                         point=r['point'], agent=r['agent'], r1=r1s, r2=r2s))
    if problems:
        raise SheetError(f'{path}:\n  ' + '\n  '.join(problems))

    lib = []
    workdir = Path(workdir)
    for (series, rep), g in pd.DataFrame(rows).groupby(['series', 'rep'], sort=False):
        name = f'{series}_rep{rep}' if has_rep else series
        concs = sorted(g['conc'].unique())
        if concs[0] != 0:
            raise SheetError(f'{path}: series {name} has no input (a row with conc 0)')
        if len(concs) < 2:
            raise SheetError(f'{path}: series {name} has only an input; add rows with conc > 0')
        tit = []
        for k, (conc, h) in enumerate(g.groupby('conc', sort=True)):
            if h['fraction'].nunique() > 1:
                raise SheetError(f'{path}: series {name}, conc {conc}: rows give different fractions')
            points = set(h['point']) - {''}
            tnum = int(points.pop()) if points else k
            tit.append((tnum, h['conc_text'].iloc[0], h['frac_text'].iloc[0]))
            lib.append(dict(name=f'{name}_{tnum}', series=name, group=series, cell=cell, tnum=tnum,
                            deep=bool(nucleosomes), agent=h['agent'].iloc[0],
                            fastq1=','.join(f for fs in h['r1'] for f in fs),
                            fastq2=','.join(f for fs in h['r2'] for f in fs),
                            titration=str((workdir / 'titration' / f'{name}.tsv').resolve())))
        # as written: rounding f to 6 digits ({f:g}) shifted GEO-reproduction scores by 1e-5
        _write_if_changed(workdir / 'titration' / f'{name}.tsv',
                          'tnum\tconc\tsurvival_fraction\n' + ''.join(f'{t}\t{c}\t{f}\n' for t, c, f in sorted(tit)))
    short = sorted({r['frac_text'] for r in rows if float(r['frac_text']) < 1 and _sig_digits(r['frac_text']) <= 6})
    if note and short:
        note(f'Note: {path}: {len(short)} soluble fraction(s) have 6 or fewer significant digits (e.g. {short[0]}). '
             'They are used as written. If they were rounded (e.g. by a spreadsheet), give the full values: '
             "Park's GEO tables used full precision, and 6-digit values change scores by up to 1e-5.")
    return _finish(pd.DataFrame(lib), path)


def read_config_sheet(samples, titration_dir, resolve):
    """The GEO-style sample sheet (cell, sample, agent, tnum, rep, depth, fastq1, fastq2)."""
    S = pd.read_csv(resolve(samples), sep='\t', dtype=str, comment='#').fillna('')
    S['tnum'] = S['tnum'].astype(int)
    S['deep'] = S['depth'].str.lower().eq('deep')
    deep = S['deep'].map({True: '_deep', False: ''})
    csa = S['cell'] + '_' + S['sample'] + '_' + S['agent']
    S['series'] = csa + '_' + S['rep'] + 'rep' + deep
    S['name'] = csa + '_' + S['tnum'].astype(str) + '_' + S['rep'] + 'rep' + deep
    S['group'] = csa + deep
    S['titration'] = [resolve(Path(titration_dir) / f'{x}.tsv') for x in csa]
    for col in ('fastq1', 'fastq2'):
        S[col] = [','.join(resolve(f) for f in v.split(',')) for v in S[col]]
    return _finish(S, resolve(samples))


def _finish(S, path):
    if not S['name'].is_unique:
        dup = sorted(S.loc[S['name'].duplicated(), 'name'])
        raise SheetError(f'{path}: duplicate libraries {dup}')
    return S[['name', 'series', 'group', 'cell', 'tnum', 'deep', 'agent', 'fastq1', 'fastq2',
              'titration']].set_index('name')


def _write_if_changed(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


def find_fastqs(dirs):
    """Read-1 FASTQ files in the given folders, for `condense init`."""
    out = []
    for d in dirs:
        for f in sorted(Path(d).expanduser().iterdir()):
            if re.search(r'\.f(ast)?q(\.gz)?$', f.name) and any(re.search(p, f.name) for p, _ in R1_PATTERNS):
                out.append(f)
    return out

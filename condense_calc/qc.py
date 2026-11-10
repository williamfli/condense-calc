"""Quality control: per-sample read statistics, per-series checks, replicate agreement, and a report.

Levels (config `qc` or `condense-calc run --qc`):
  off      no QC
  default  alignment and read-filter summary, fragment lengths, titration, score distributions,
           score vs input/supernatant counts, c-half fit success, nucleosome spacing (deep),
           replicate agreement; flags; one figure per series and per replicate group
  verbose  adds duplicate rate, GC dependence, a chromosome-wide score profile and the
           dinucleotide periodicity around called dyads

Several checks follow Park et al. 2025 (fragment lengths: ED Fig. 2c; input vs supernatant:
ED Fig. 2d,f; dinucleotide periodicity: ED Fig. 2e; replicate correlation: ED Fig. 2j).
"""
import html
import json
import math
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from . import titration

TLEN_MAX = 1000
MONO = (126, 184)            # the peak caller's mononucleosome length window
NAME_RE = re.compile(r'_(\d+)_(\d+)rep(_deep)?$')

THRESHOLDS = {
    'min_alignment_rate': (0.80, 'bowtie2 overall alignment rate'),
    'min_used_fraction': (0.50, 'share of read pairs kept as fragments (proper pair, edit distance, MAPQ >= min_mapq)'),
    'max_lowmapq_fraction': (0.15, 'share of fragments with MAPQ < 10 (multi-mapping), before the MAPQ filter'),
    'min_mono_fraction': (0.50, f'share of fragments {MONO[0]}-{MONO[1]} bp'),
    'min_tlen_mode': (130, 'fragment-length mode, lower bound'),
    'max_tlen_mode': (170, 'fragment-length mode, upper bound'),
    'max_dup_fraction': (0.30, 'share of duplicate fragments (verbose only)'),
    'max_zero_input_fraction': (0.15, 'share of bins without input reads'),
    'max_score_input_rho': (0.30, '|Spearman| of the binned score with the input read count'),
    'min_replicate_rho': (0.60, 'Spearman between replicates at a titration point'),
    'min_fit_fraction': (0.90, 'share of bins with a c-half fit (R^2 above the cutoff)'),
}


def tnum_of(name):
    """Titration point of a library name: GEO style (…_<tnum>_<rep>rep[_deep]) or sheet style (…_<point>)."""
    m = NAME_RE.search(name) or re.search(r'_(\d+)$', name)
    return int(m.group(1)) if m else None


def _write_json(path, d):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as fh:
        json.dump(d, fh, indent=1, default=float)


def _spearman(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float('nan')
    return float(pd.Series(a[ok]).rank().corr(pd.Series(b[ok]).rank()))


# ---------------------------------------------------------------- per sample

def parse_bowtie2_log(path):
    txt = Path(path).read_text()

    def num(pat, cast=int):
        m = re.search(pat, txt)
        return cast(m.group(1)) if m else None
    return {
        'read_pairs': num(r'(\d+) reads; of these'),
        'concordant_once': num(r'(\d+) \([\d.]+%\) aligned concordantly exactly 1 time'),
        'concordant_multi': num(r'(\d+) \([\d.]+%\) aligned concordantly >1 times'),
        'alignment_rate': (num(r'([\d.]+)% overall alignment rate', float) or 0) / 100,
    }


def sample_stats(bam, bowtie2_log, mm=10, min_mapq=10, threads=1, verbose=False):
    """Reads the fragments the pipeline counts (forward read of a proper pair, edit distance
    <= mm, MAPQ >= min_mapq) once with samtools, and summarizes them. lowmapq_fraction is taken
    before the MAPQ filter, so it shows how much of the library maps to several places."""
    d = {'sample': Path(bam).name.split('.')[0], **parse_bowtie2_log(bowtie2_log)}
    cmd = (f"set -o pipefail; samtools view -@ {max(threads - 1, 1)} -f 0x2 -F 0x1C -e '[NM]<={mm}' "
           f"'{bam}' | cut -f 3,4,5,9")
    proc = subprocess.Popen(['bash', '-c', cmd], stdout=subprocess.PIPE)
    hist = np.zeros(TLEN_MAX + 2, dtype=np.int64)       # last bin: longer than TLEN_MAX
    n = n_all = lowq = dup = 0

    def add(df):
        nonlocal n, n_all, lowq, dup
        n_all += len(df)
        lowq += int((df['mapq'] < 10).sum())
        df = df[df['mapq'] >= min_mapq]
        n += len(df)
        hist[:] += np.bincount(np.minimum(df['tlen'].abs().to_numpy(), TLEN_MAX + 1), minlength=TLEN_MAX + 2)
        if verbose:
            dup += int(df.duplicated(['chrom', 'pos', 'tlen']).sum())

    reader = pd.read_csv(proc.stdout, sep='\t', header=None, names=['chrom', 'pos', 'mapq', 'tlen'],
                         dtype={'chrom': str, 'pos': np.int64, 'mapq': np.int64, 'tlen': np.int64},
                         chunksize=5_000_000)
    carry = None
    for chunk in reader:       # the BAM is sorted; keep the last position for the next chunk
        if carry is not None:
            chunk = pd.concat([carry, chunk], ignore_index=True)
        tail = (chunk['chrom'] == chunk['chrom'].iat[-1]) & (chunk['pos'] == chunk['pos'].iat[-1])
        carry = chunk[tail]
        add(chunk[~tail])
    if carry is not None:
        add(carry)
    if proc.wait() != 0:
        raise RuntimeError(f'samtools failed on {bam}')
    tl = np.arange(TLEN_MAX + 2)
    frac = lambda m: float(hist[m].sum() / max(n, 1))
    cum = np.cumsum(hist)
    d.update({
        'fragments': n,
        'used_fraction': n / d['read_pairs'] if d['read_pairs'] else float('nan'),
        'min_mapq': min_mapq,
        'lowmapq_fragments': lowq,
        'lowmapq_fraction': lowq / max(n_all, 1),
        'tlen_median': int(np.searchsorted(cum, n / 2)) if n else None,
        'tlen_mode': int(50 + np.argmax(hist[50:TLEN_MAX + 1])) if n else None,
        'mono_fraction': frac((tl >= MONO[0]) & (tl <= MONO[1])),
        'subnucleosomal_fraction': frac(tl < 120),
        'long_fraction': frac(tl > 250),
        'tlen_hist': hist[:TLEN_MAX + 1].tolist(),
    })
    if verbose:
        d['dup_fraction'] = dup / max(n, 1)
    return d


# ---------------------------------------------------------------- per series

def _read_table(path, usecols=None):
    return pd.read_csv(path, sep='\t', usecols=usecols, na_values=['NA', 'nan'])


def series_stats(series, sample_jsons, bins, table, titration_path, fit=None, peaks=(), fasta=None,
                 verbose=False):
    samples = [json.loads(Path(p).read_text()) for p in sample_jsons]
    tit = titration.read(titration_path)
    tn = sorted(tit.frac)
    d = {'series': series, 'samples': [s['sample'] for s in samples],
         'titration': {'tnum': tn, 'conc': [tit.conc[t] for t in tn], 'frac': [tit.frac[t] for t in tn]},
         'titration_monotonic': bool(all(np.diff([tit.frac[t] for t in tn]) <= 1e-9))}

    # bin counts of every sample, joined on bin
    key = ['Chromosome', 'Start']
    counts, gc = None, None
    for p in bins:
        t = _read_table(p)
        name = t.columns[3]
        if gc is None and 'GCcontent' in t.columns:
            gc = t[key + ['GCcontent']]
        t = t[key + [name]]
        counts = t if counts is None else counts.merge(t, on=key, how='outer')
    inp = next(c for c in counts.columns[2:] if tnum_of(c) == 0)
    c0 = counts[inp].fillna(0).to_numpy(float)
    d['bins'] = len(counts)
    d['zero_input_fraction'] = float((c0 == 0).mean())
    d['input_cv'] = float(np.std(c0[c0 > 0]) / np.mean(c0[c0 > 0])) if (c0 > 0).any() else float('nan')

    sc = _read_table(table)
    m = sc.merge(counts, on=key, how='left', suffixes=('', '_count'))
    if gc is not None:
        m = m.merge(gc, on=key, how='left')
    pts = []
    for col in sc.columns[3:]:
        s = m[col].to_numpy(float)
        cnt_col = col + '_count' if col + '_count' in m.columns else None
        row = {'name': col, 'tnum': tnum_of(col), 'bins_scored': int(np.isfinite(s).sum()),
               'median': float(np.nanmedian(s)), 'q05': float(np.nanquantile(s, 0.05)),
               'q95': float(np.nanquantile(s, 0.95)),
               'rho_input': _spearman(s, m[inp].to_numpy(float)),
               'rho_supernatant': _spearman(s, m[cnt_col].to_numpy(float)) if cnt_col else float('nan')}
        if verbose and gc is not None:
            row['rho_gc'] = _spearman(s, m['GCcontent'].to_numpy(float))
        pts.append(row)
    d['points'] = pts

    if fit:
        f = _read_table(fit)
        ch = f['C-half'].to_numpy(float)
        d['fit'] = {'bins': len(f), 'fitted_fraction': float(np.isfinite(ch).mean()),
                    'chalf_median': float(np.nanmedian(ch)) if np.isfinite(ch).any() else float('nan'),
                    'chalf_hist': np.histogram(ch[np.isfinite(ch)], bins=50)[0].tolist(),
                    'chalf_edges': np.histogram(ch[np.isfinite(ch)], bins=50)[1].tolist(),
                    'rsq_hist': np.histogram(f['R-squared'].to_numpy(float), bins=50, range=(0, 1))[0].tolist()}
    if peaks:
        n, spacing = 0, np.zeros(401, dtype=np.int64)
        for p in peaks:
            pos = np.sort(_read_table(p, usecols=[1]).iloc[:, 0].to_numpy(np.int64))
            n += len(pos)
            spacing += np.bincount(np.minimum(np.diff(pos), 400), minlength=401)
        sm = np.convolve(spacing[:400], np.ones(9) / 9, mode='same')
        d['nucleosomes'] = {'called': int(n), 'spacing_hist': spacing.tolist(),
                            'spacing_median': int(np.searchsorted(np.cumsum(spacing), spacing.sum() / 2)),
                            'spacing_mode': int(140 + np.argmax(sm[140:260]))}
    if verbose:
        d['profile'] = _chrom_profile(m, sc.columns[3:])
        if peaks and fasta:
            d['periodicity'] = _dinucleotide_periodicity(peaks[0], fasta)
    return d


def _chrom_profile(m, cols, bin_size=1_000_000):
    chrom = m['Chromosome'].iat[0]
    x = m[m['Chromosome'] == chrom]
    g = x.groupby(x['Start'] // bin_size)[list(cols)].mean()
    return {'chrom': chrom, 'mb': g.index.tolist(), 'scores': {c: g[c].tolist() for c in cols}}


def _dinucleotide_periodicity(peak_file, fasta, half=80, max_dyads=200_000):
    """AA/AT/TA/TT and CC/CG/GC/GG frequency at each position around dyads (ED Fig. 2e)."""
    import pysam
    t = _read_table(peak_file)
    if len(t) > max_dyads:
        t = t.sample(max_dyads, random_state=0)
    ww, ss = np.zeros(2 * half), np.zeros(2 * half)
    k = 0
    with pysam.FastaFile(fasta) as fa:
        for chrom, pos in zip(t.iloc[:, 0], t.iloc[:, 1]):
            if pos - half < 0:
                continue
            seq = fa.fetch(chrom, pos - half, pos + half + 1).upper()
            if len(seq) != 2 * half + 1:
                continue
            di = np.array([seq[i:i + 2] for i in range(2 * half)])
            ww += np.isin(di, ['AA', 'AT', 'TA', 'TT'])
            ss += np.isin(di, ['CC', 'CG', 'GC', 'GG'])
            k += 1
    return {'dyads': k, 'offset': list(range(-half, half)), 'ww': (ww / max(k, 1)).tolist(),
            'ss': (ss / max(k, 1)).tolist()}


# ---------------------------------------------------------------- replicates

def replicate_stats(group, tables_by_kb):
    """Spearman and Pearson between replicate series at each shared titration point."""
    d = {'group': group, 'pairs': []}
    for kb, paths in tables_by_kb.items():
        tabs = [_read_table(p) for p in paths]
        for i in range(len(tabs)):
            for j in range(i + 1, len(tabs)):
                a, b = tabs[i], tabs[j]
                ca = {tnum_of(c): c for c in a.columns[3:]}
                cb = {tnum_of(c): c for c in b.columns[3:]}
                m = a.merge(b, on=['Chromosome', 'Start'], suffixes=('_a', '_b'))
                for t in sorted(set(ca) & set(cb)):
                    xa = m[ca[t] + ('_a' if ca[t] == cb[t] else '')].to_numpy(float)
                    xb = m[cb[t] + ('_b' if ca[t] == cb[t] else '')].to_numpy(float)
                    ok = np.isfinite(xa) & np.isfinite(xb)
                    row = {'kb': kb, 'a': Path(paths[i]).name, 'b': Path(paths[j]).name, 'tnum': t,
                           'bins': int(ok.sum()), 'spearman': _spearman(xa, xb),
                           'pearson': float(np.corrcoef(xa[ok], xb[ok])[0, 1]) if ok.sum() > 2 else float('nan')}
                    if t == max(set(ca) & set(cb)):
                        rng = np.random.default_rng(0)
                        idx = rng.choice(np.flatnonzero(ok), size=min(50_000, ok.sum()), replace=False)
                        row['sample_xy'] = [xa[idx].round(4).tolist(), xb[idx].round(4).tolist()]
                    d['pairs'].append(row)
    return d


# ---------------------------------------------------------------- figures

def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def series_figure(d, samples, png):
    plt = _plt()
    panels = ['reads', 'tlen', 'titration', 'scores', 'rho']
    if 'fit' in d:
        panels += ['chalf', 'rsq']
    if 'nucleosomes' in d:
        panels.append('spacing')
    if 'profile' in d:
        panels.append('profile')
    if any('rho_gc' in p for p in d['points']):
        panels.append('gc')
    if 'periodicity' in d:
        panels.append('periodicity')
    ncol = 3
    nrow = math.ceil(len(panels) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.6 * nrow), squeeze=False)
    tn = [tnum_of(s['sample']) for s in samples]
    cmap = plt.get_cmap('viridis')
    color = {t: cmap(i / max(len(tn) - 1, 1)) for i, t in enumerate(sorted(tn))}
    for ax, panel in zip(axes.flat, panels):
        if panel == 'reads':
            x = np.arange(len(samples))
            ax.bar(x, [s['read_pairs'] / 1e6 for s in samples], color='0.8', label='read pairs')
            ax.bar(x, [s['fragments'] / 1e6 for s in samples], color='C0', label='fragments kept')
            ax.bar(x, [s.get('lowmapq_fragments', s['fragments'] * s['lowmapq_fraction']) / 1e6 for s in samples],
                   color='C3', label='MAPQ < 10' + (' (not kept)' if samples[0].get('min_mapq', 0) >= 10 else ''))
            ax.set_xticks(x, [f"#{t}" for t in tn])
            ax.set_ylabel('millions')
            ax.set_title('Reads per titration point')
            ax.legend(fontsize=7)
        elif panel == 'tlen':
            for s, t in zip(samples, tn):
                h = np.asarray(s['tlen_hist'][:401], float)
                ax.plot(np.arange(401), h / max(h.sum(), 1), color=color[t], lw=1, label=f'#{t}')
            ax.axvspan(*MONO, color='0.9', zorder=0)
            ax.set_xlabel('fragment length (bp)')
            ax.set_title('Fragment lengths')
            ax.legend(fontsize=6, ncol=2)
        elif panel == 'titration':
            t = d['titration']
            ax.plot(t['conc'], t['frac'], 'o-', color='0.4')
            for c, f, k in zip(t['conc'], t['frac'], t['tnum']):
                if k in tn:
                    ax.annotate(f'#{k}', (c, f), fontsize=7, xytext=(3, 3), textcoords='offset points')
            ax.set_xlabel('concentration')
            ax.set_ylabel('soluble fraction f')
            ax.set_title('Titration (input values)')
        elif panel == 'scores':
            p = d['points']
            ax.errorbar([q['tnum'] for q in p], [q['median'] for q in p],
                        yerr=[[q['median'] - q['q05'] for q in p], [q['q95'] - q['median'] for q in p]],
                        fmt='o', color='C0')
            ax.set_xlabel('titration point')
            ax.set_ylabel('score (median, 5-95%)')
            ax.set_title('Score distributions')
        elif panel == 'rho':
            p = d['points']
            ax.plot([q['tnum'] for q in p], [q['rho_supernatant'] for q in p], 'o-', label='vs supernatant count')
            ax.plot([q['tnum'] for q in p], [q['rho_input'] for q in p], 's-', label='vs input count')
            ax.axhline(0, color='0.7', lw=0.8)
            ax.set_ylim(-1, 1)
            ax.set_xlabel('titration point')
            ax.set_ylabel('Spearman')
            ax.set_title('Score vs read counts')
            ax.legend(fontsize=7)
        elif panel == 'chalf':
            f = d['fit']
            e = np.asarray(f['chalf_edges'])
            ax.bar(e[:-1], f['chalf_hist'], width=np.diff(e), align='edge', color='C2')
            ax.set_xlabel('c-half')
            ax.set_title(f"c-half ({f['fitted_fraction']:.1%} of bins fitted)")
        elif panel == 'rsq':
            ax.bar(np.linspace(0, 1, 51)[:-1], d['fit']['rsq_hist'], width=0.02, align='edge', color='C2')
            ax.set_xlabel('R^2')
            ax.set_title('c-half fit quality')
        elif panel == 'spacing':
            nd = d['nucleosomes']
            ax.plot(np.arange(400), nd['spacing_hist'][:400], color='C4')
            ax.set_xlabel('distance to next dyad (bp)')
            ax.set_title(f"{nd['called']:,} nucleosomes; spacing mode {nd['spacing_mode']} bp")
        elif panel == 'profile':
            pr = d['profile']
            for c, v in pr['scores'].items():
                ax.plot(pr['mb'], v, '.-', color=color.get(tnum_of(c), 'k'), lw=0.8, ms=3, label=f'#{tnum_of(c)}')
            ax.set_xlabel(f"{pr['chrom']} (Mb)")
            ax.set_ylabel('mean score, 1 Mb')
            ax.set_title('Chromosome profile')
        elif panel == 'gc':
            p = d['points']
            ax.plot([q['tnum'] for q in p], [q['rho_gc'] for q in p], 'o-', color='C5')
            ax.axhline(0, color='0.7', lw=0.8)
            ax.set_xlabel('titration point')
            ax.set_ylabel('Spearman, score vs GC')
            ax.set_title('GC dependence')
        elif panel == 'periodicity':
            pe = d['periodicity']
            ax.plot(pe['offset'], pe['ww'], label='AA/AT/TA/TT')
            ax.plot(pe['offset'], pe['ss'], label='CC/CG/GC/GG')
            ax.set_xlabel('position from dyad (bp)')
            ax.set_title(f"Dinucleotides at {pe['dyads']:,} dyads")
            ax.legend(fontsize=7)
    from matplotlib.ticker import MaxNLocator
    for ax, panel in zip(axes.flat, panels):
        if panel in ('scores', 'rho', 'gc'):
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for ax in list(axes.flat)[len(panels):]:
        ax.axis('off')
    fig.suptitle(d['series'])
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)


def replicate_figure(d, png):
    plt = _plt()
    kbs = sorted({p['kb'] for p in d['pairs']}, key=lambda k: int(k[:-2]))
    fig, axes = plt.subplots(len(kbs), 2, figsize=(10, 3.8 * len(kbs)), squeeze=False)
    for (a1, a2), kb in zip(axes, kbs):
        rows = [p for p in d['pairs'] if p['kb'] == kb]
        for pair in sorted({(p['a'], p['b']) for p in rows}):
            r = [p for p in rows if (p['a'], p['b']) == pair]
            a1.plot([p['tnum'] for p in r], [p['spearman'] for p in r], 'o-', label=' vs '.join(pair))
        a1.set_ylim(0, 1)
        a1.set_xlabel('titration point')
        a1.set_ylabel('Spearman')
        a1.set_title(f'Replicate agreement, {kb}')
        a1.legend(fontsize=6)
        xy = next((p for p in rows if 'sample_xy' in p), None)
        if xy:
            a2.hexbin(*xy['sample_xy'], gridsize=60, bins='log', mincnt=1)
            a2.set_xlabel(xy['a'])
            a2.set_ylabel(xy['b'])
            a2.set_title(f"point #{xy['tnum']}: Spearman {xy['spearman']:.2f}")
    fig.suptitle(d['group'])
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------- report

def flags(samples, series, reps, th):
    out = []

    def flag(level, name, metric, value, limit):
        out.append({'level': level, 'name': name, 'metric': metric, 'value': value, 'threshold': limit,
                    'meaning': THRESHOLDS.get(limit.split()[0], ('', ''))[1] if limit else ''})
    for s in samples:
        n = s['sample']
        checks = [('alignment_rate', '<', 'min_alignment_rate'), ('used_fraction', '<', 'min_used_fraction'),
                  ('lowmapq_fraction', '>', 'max_lowmapq_fraction'), ('mono_fraction', '<', 'min_mono_fraction'),
                  ('tlen_mode', '<', 'min_tlen_mode'), ('tlen_mode', '>', 'max_tlen_mode'),
                  ('dup_fraction', '>', 'max_dup_fraction')]
        for metric, op, key in checks:
            v = s.get(metric)
            if v is not None and (v < th[key] if op == '<' else v > th[key]):
                flag('warn', n, metric, v, f'{key} = {th[key]}')
    for d in series:
        if not d['titration_monotonic']:
            flag('warn', d['series'], 'titration', 'soluble fraction rises with tnum', '')
        if d['zero_input_fraction'] > th['max_zero_input_fraction']:
            flag('warn', d['series'], 'zero_input_fraction', d['zero_input_fraction'],
                 f"max_zero_input_fraction = {th['max_zero_input_fraction']}")
        frac = dict(zip(d['titration']['tnum'], d['titration']['frac']))
        for p in d['points']:
            if abs(p['rho_input']) > th['max_score_input_rho']:
                # where most nucleosomes are still soluble the score is mostly counting noise
                flag('info' if frac.get(p['tnum'], 0) > 0.5 else 'warn', p['name'], 'rho_input', p['rho_input'],
                     f"max_score_input_rho = {th['max_score_input_rho']}")
        if 'fit' in d and d['fit']['fitted_fraction'] < th['min_fit_fraction']:
            flag('warn', d['series'], 'fitted_fraction', d['fit']['fitted_fraction'],
                 f"min_fit_fraction = {th['min_fit_fraction']}")
    for d in reps:
        for p in d['pairs']:
            if p['spearman'] < th['min_replicate_rho']:
                flag('info' if p['tnum'] <= 2 else 'warn', f"{d['group']} #{p['tnum']} {p['kb']}",
                     'replicate_spearman', p['spearman'], f"min_replicate_rho = {th['min_replicate_rho']}")
    return out


def _fmt(v):
    if isinstance(v, float):
        return f'{v:.3g}' if abs(v) < 1e4 else f'{v:,.0f}'
    if isinstance(v, int):
        return f'{v:,}'
    return str(v)


def report(sample_jsons, series_jsons, rep_jsons, overrides, outdir, ed2_jsons=()):
    outdir = Path(outdir)
    load = lambda ps: [json.loads(Path(p).read_text()) for p in ps]
    samples, series, reps = load(sample_jsons), load(series_jsons), load(rep_jsons)
    ed2 = load(ed2_jsons)
    th = {k: v[0] for k, v in THRESHOLDS.items()}
    for kv in overrides:
        k, v = kv.split('=', 1)
        if k not in th:
            raise ValueError(f'unknown QC threshold {k!r}; known: {", ".join(th)}')
        th[k] = float(v)
    fl = flags(samples, series, reps, th)

    sample_cols = ['read_pairs', 'alignment_rate', 'fragments', 'used_fraction', 'lowmapq_fraction',
                   'tlen_mode', 'tlen_median', 'mono_fraction', 'subnucleosomal_fraction', 'dup_fraction']
    st = pd.DataFrame([{'sample': s['sample'], **{c: s.get(c) for c in sample_cols}} for s in samples])
    st = st.dropna(axis=1, how='all')
    rows = []
    for d in series:
        base = {'series': d['series'], 'bins': d['bins'], 'zero_input_fraction': d['zero_input_fraction'],
                'input_cv': d['input_cv']}
        if 'fit' in d:
            base.update(fitted_fraction=d['fit']['fitted_fraction'], chalf_median=d['fit']['chalf_median'])
        if 'nucleosomes' in d:
            base.update(nucleosomes=d['nucleosomes']['called'], spacing_mode=d['nucleosomes']['spacing_mode'])
        for p in d['points']:
            rows.append({**base, 'point': p['name'], 'tnum': p['tnum'], 'bins_scored': p['bins_scored'],
                         'score_median': p['median'], 'rho_input': p['rho_input'],
                         'rho_supernatant': p['rho_supernatant'], **({'rho_gc': p['rho_gc']} if 'rho_gc' in p else {})})
    pt = pd.DataFrame(rows)
    rt = pd.DataFrame([{k: v for k, v in p.items() if k != 'sample_xy'} | {'group': d['group']}
                       for d in reps for p in d['pairs']])

    long = []
    for df, level, key in [(st, 'sample', 'sample'), (pt, 'point', 'point')]:
        if len(df):
            for _, r in df.iterrows():
                long += [{'level': level, 'name': r[key], 'metric': c, 'value': r[c]}
                         for c in df.columns if c not in (key,) and pd.notna(r[c])]
    for _, r in rt.iterrows():
        long.append({'level': 'replicates', 'name': f"{r['group']} {r['kb']} #{r['tnum']}",
                     'metric': 'spearman', 'value': r['spearman']})
    pd.DataFrame(long).to_csv(outdir / 'summary.tsv', sep='\t', index=False)
    pd.DataFrame(fl, columns=['level', 'name', 'metric', 'value', 'threshold', 'meaning']).to_csv(
        outdir / 'flags.tsv', sep='\t', index=False)

    def table(df):
        if not len(df):
            return '<p>none</p>'
        head = ''.join(f'<th>{html.escape(str(c))}</th>' for c in df.columns)
        body = ''.join('<tr>' + ''.join(f'<td>{html.escape(_fmt(v))}</td>' for v in r) + '</tr>'
                       for r in df.itertuples(index=False))
        return f'<table><tr>{head}</tr>{body}</table>'
    parts = ['<!doctype html><meta charset="utf-8"><title>condense-calc QC</title>',
             '<style>body{font:14px sans-serif;margin:24px;max-width:1500px}table{border-collapse:collapse;'
             'margin:8px 0 20px}td,th{border:1px solid #ccc;padding:3px 7px;text-align:right}'
             'th{background:#f3f3f3}img{max-width:100%;border:1px solid #ddd}</style>',
             '<h1>condense-calc QC report</h1>',
             f'<h2>Flags ({len(fl)})</h2>',
             table(pd.DataFrame(fl)) if fl else '<p>No flags.</p>',
             '<h2>Samples</h2>', table(st),
             '<h2>Titration points</h2>', table(pt.drop(columns=[c for c in ('bins',) if c in pt]))]
    for d in series:
        parts += [f"<h2>{html.escape(d['series'])}</h2>", f"<img src='series/{html.escape(d['series'])}.png'>"]
    if reps:
        parts += ['<h2>Replicates</h2>', table(rt)]
        parts += [f"<img src='replicates/{html.escape(d['group'])}.png'>" for d in reps]
    if ed2:
        rows = [{'series': e['series'], 'chrom': e['chrom'], 'panels': ''.join(e['panels']),
                 'd_mask': e['d']['mask'],
                 'e_ww_period_bp': e.get('e', {}).get('ww_period_bp'),
                 'f_rho_score_input': e.get('f', {}).get('spearman_score_input'),
                 'f_rho_score_sup': e.get('f', {}).get('spearman_score_sup'),
                 'states': len(e['states']['used']) if 'states' in e else None,
                 'h_max_cohen_d': e.get('h', {}).get('cohen_d_max')} for e in ed2]
        parts += ['<h2>Panels after Extended Data Fig. 2 (Park et al. 2025)</h2>',
                  '<p>c fragment lengths; d NCP number along the chromosome (25 kb, smoothed) with the condensability '
                  'of the last point; deep series: e dinucleotides around the dyads, f input vs last point per '
                  'nucleosome; with chromatin states: g NCP number, h called nucleosomes vs sliding windows, '
                  "i Cohen's d and Welch p between states.</p>",
                  table(pd.DataFrame(rows).dropna(axis=1, how='all'))]
        for e in ed2:
            d = f"ed2/{html.escape(e['series'])}"
            parts += [f"<h3>{html.escape(e['series'])}</h3>",
                      ''.join(f"<img src='{d}/ed2{k}.png' style='max-height:{360 if k != 'd' else 300}px;"
                              f"{'' if k != 'd' else 'display:block;'}margin:4px'>" for k in e['panels'])]
    parts += ['<h2>Thresholds</h2>',
              table(pd.DataFrame([{'name': k, 'value': th[k], 'meaning': v[1]} for k, v in THRESHOLDS.items()]))]
    (outdir / 'report.html').write_text('\n'.join(parts))
    return fl

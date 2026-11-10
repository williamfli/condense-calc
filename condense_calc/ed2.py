"""Single-dataset QC panels after Extended Data Fig. 2 of Park et al. 2025, computed and drawn as in their
notebooks (condense-seq repository, ipython_notebooks/[Extended Data Fig.2*]), from one titration series:

  c  fragment-length distribution of each sample                                          every series
  d  NCP number along one chromosome (25 kb bins, 81-bin moving mean, log10, centred),
     with the condensability of the last point                                           every series
  e  AA/AT/TA/TT and CC/CG/GC/GG frequency in 147 bp around the called dyads                deep series
  f  per nucleosome: NCP number of input vs last point, coloured by condensability          deep series
  g  NCP number of input and last point per chromatin state                                 deep, with states
  h  condensability per state: called nucleosomes vs 171 bp sliding windows                 deep, with states
  i  Cohen's d and Welch's t-test p-value between states                                    deep, with states

Park's panels come from H1 (c, d: the shallow series, ten points; e-i: deep replicate 1, chromosome 1). Here each
series gets the panels its data allow, on chromosome 1 (or the first chromosome). Differences from the notebooks:
when no cytogenetic bands are given, panel d masks assembly gaps (25 kb bins more than half N) instead of the acen
and gvar bands; a series' bins are summed to 10 kb for panel d when its bins are smaller (the same as the
pipeline's 10 kb bins); hand-drawn marks in the published figure (d, e titles, f arrows, h n.s. marks) are drawn
from the computed values.
"""
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from . import annotations, titration
from .qc import tnum_of, _spearman

MNUM = 1.6e12                    # NCPnum_ver2.py: input molecules
UNITS = {'sp': 'mM', 'spd': 'mM', 'CoH': 'mM', 'Mg': 'mM', 'Ca': 'mM'}
KEY = ['Chromosome', 'Start']


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return matplotlib, plt


def _read(path, usecols=None):
    return pd.read_csv(path, sep='\t', usecols=usecols, na_values=['NA', 'nan'])


def conc_label(agent, conc):
    return (f'[{agent}]=%.2f{UNITS.get(agent, "")}' if agent else 'conc=%.2f') % conc


# ---------------------------------------------------------------- statistics

def cohen_abs(x, y):
    """|Cohen's d| with the pooled standard deviation (statis_edit.get_CohenD; np.std, ddof 0)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    s = math.sqrt(((len(x) - 1) * x.std() ** 2 + (len(y) - 1) * y.std() ** 2) / (len(x) + len(y) - 2))
    return abs(x.mean() - y.mean()) / s if s > 0 else float('nan')


def welch_p(x, y):
    """Two-sided Welch's t-test p-value (scipy.stats.ttest_ind, equal_var=False)."""
    from scipy import stats
    if len(x) < 2 or len(y) < 2:
        return float('nan')
    return float(stats.ttest_ind(x, y, equal_var=False)[1])


def categorize(pos, states):
    """State containing each position, half-open [start, end) (Interval_dict.double_hash.find); None outside."""
    st = np.array([s for s, _, _ in states], dtype=np.int64)
    ed = np.array([e for _, e, _ in states], dtype=np.int64)
    name = np.array([n for _, _, n in states], dtype=object)
    pos = np.asarray(pos, dtype=np.int64)
    if not len(st):
        return np.full(len(pos), None, dtype=object)
    i = np.searchsorted(st, pos, side='right') - 1
    j = np.clip(i, 0, None)
    ok = (i >= 0) & (pos < ed[j])
    return np.where(ok, name[j], None)


# ---------------------------------------------------------------- panel d helpers

def rebin25(v, start, length, step=25000):
    """Overlap-weighted mean of 10 kb values over 25 kb bins k*25 kb (rbin_data_mean): each 10 kb value spread on
    a 5 kb grid; each 25 kb bin the mean of its five grid cells that have a value."""
    n5 = (length // step + 1) * 5
    g = np.full(n5, np.nan)
    i = (np.asarray(start) // 5000).astype(int)
    for k in (0, 1):
        j = i + k
        m = j < n5
        g[j[m]] = np.asarray(v, float)[m]
    with warnings.catch_warnings():                          # all-NaN bins stay NaN
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(g.reshape(-1, 5), axis=1)


def smooth(x, win=81):
    """Centred moving mean over win values that skips NaN, truncated at the ends (slow_moving_average2); also for
    series shorter than win."""
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    cs = np.concatenate([[0.0], np.cumsum(np.where(ok, x, 0.0))])
    cn = np.concatenate([[0], np.cumsum(ok)])
    i = np.arange(len(x))
    lo, hi = np.clip(i - win // 2, 0, len(x)), np.clip(i + win // 2 + 1, 0, len(x))
    s, n = cs[hi] - cs[lo], cn[hi] - cn[lo]
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(n > 0, s / n, np.nan)


def mask_spans(raw, gap=10):
    """Merge [start, end) bin spans closer than `gap` bins (Park's d notebook)."""
    spans = []
    for a, b in sorted(raw):
        if spans and a - spans[-1][1] <= gap:
            spans[-1][1] = max(spans[-1][1], b)
        else:
            spans.append([a, b])
    return spans


def bin_score(t, c, f):
    """Park's bin score of summed counts (NCPscore_ver4.py, bins without input dropped): -ln[((t+1)/T)/((c+1)/C)] -
    ln f, totals over the bins with input."""
    ok = c > 0
    T, C = (t[ok] + 1).sum(), (c[ok] + 1).sum()
    out = np.full(len(t), np.nan)
    out[ok] = -np.log(((t[ok] + 1) / T) / ((c[ok] + 1) / C)) - np.log(f)
    return out


# ---------------------------------------------------------------- panels

def panel_c(samples, tit, agent, out, R, plt, mpl):
    fig, ax = plt.subplots(figsize=(5, 4))
    cols = mpl.colormaps['copper'](np.linspace(0.01, 0.7, 10))
    alphas = np.linspace(0.3, 1, 10)[::-1]
    R['c'] = {}
    for s in sorted(samples, key=lambda s: tnum_of(s['sample'])):
        t = tnum_of(s['sample'])
        h = np.asarray(s['tlen_hist'], float)
        k = min(t, 9)
        ax.plot(np.arange(len(h)), h, color=cols[k], alpha=alphas[k], label=conc_label(agent, tit.conc[t]))
        R['c'][t] = dict(mode=int(np.argmax(h[50:]) + 50) if h.sum() else None, fragments=int(h.sum()),
                         mono_126_184=float(h[126:185].sum() / max(h.sum(), 1)))
    ax.grid(True)
    ax.set_xticks(range(0, 500, 50), [str(x) for x in range(0, 500, 50)])
    ax.set_xlim(-25, 525)
    ax.set_xlabel('Read length (bp)'); ax.set_ylabel('Read counts'); ax.set_title('Read length distribution')
    ax.legend(loc='best', fontsize=7)
    fig.tight_layout(); fig.savefig(out / 'ed2c.png', dpi=200); plt.close(fig)


def bins_10kb(bins, bin_size):
    """{tnum: DataFrame Chromosome, Start, n} at 10 kb, summed from smaller bins when needed."""
    out = {}
    for p in bins:
        b = _read(p)
        name = b.columns[3]
        b = b[['Chromosome', 'Start', name]].rename(columns={name: 'n'})
        if bin_size != 10000:
            b = b.assign(Start=b.Start // 10000 * 10000).groupby(KEY, sort=False, as_index=False)['n'].sum()
        out[tnum_of(name)] = b
    return out


def panel_d(series, b10, tit, agent, mscale, chrom, seq, bands, sup, out, R, plt, mpl):
    L = len(seq)
    nb = L // 25000 + 1
    tracks = {}
    for t, b in sorted(b10.items()):
        num = np.round(MNUM * mscale * tit.frac[t] * b['n'] / max(b['n'].sum(), 1))
        c1 = (b.Chromosome == chrom).to_numpy()
        sig = smooth(rebin25(num[c1].to_numpy(float), b.loc[c1, 'Start'].to_numpy(), L))
        x = np.log10(sig + 1)
        tracks[t] = x - np.nanmean(x)
    c, s = b10[0], b10[sup]
    m = c.merge(s, on=KEY, how='left', suffixes=('_c', '_t')).fillna({'n_t': 0})
    sc = bin_score(m['n_t'].to_numpy(float), m['n_c'].to_numpy(float), tit.frac[sup])
    c1 = (m.Chromosome == chrom).to_numpy()
    sct = smooth(rebin25(sc[c1], m.loc[c1, 'Start'].to_numpy(), L))
    if bands:
        raw = [[a // 25000, b // 25000 + 1] for a, b in bands]
        how = 'acen and gvar bands'
    else:
        nfrac = np.array([(seq[k * 25000:(k + 1) * 25000] == ord('N')).mean() if k * 25000 < L else 1.0
                          for k in range(nb)])
        raw = [[k, k + 1] for k in np.flatnonzero(nfrac > 0.5)]
        how = 'assembly gaps'
    spans = mask_spans(raw)
    fig, ax = plt.subplots(figsize=(15, 3))
    ax2 = ax.twinx()
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    jet = mpl.colormaps['jet_r'](np.linspace(0.1, 0.9, 10))
    xs = np.arange(nb)
    for t in tracks:
        ax.plot(xs, tracks[t], color=jet[min(t, 9)], label=conc_label(agent, tit.conc[t]))
    ax2.plot(xs, sct, color='k', alpha=0.5, linestyle='--')
    for a, b in spans:
        ax.axvspan(a + 50 if (a > 0 and b - a > 51) else a, b - 1, facecolor='white', zorder=2.5, lw=0)
    step = 10 if L > 50_000_000 else max(1, L // 5_000_000) if L > 5_000_000 else None
    if step:
        ticks = [(k, str(mb)) for mb in range(0, L // 1_000_000 + 1, step)
                 if 0 <= (k := int(round(mb * 1e6 / 25000 - 0.5))) < nb]
        ax.set_xticks([k for k, _ in ticks], [lab for _, lab in ticks], fontsize=10)
        ax.set_xlabel('Position (Mb)', fontsize=15)
    else:
        ax.set_xlabel('Position (25 kb bins)', fontsize=15)
    ax.set_xlim(0, nb + 1)
    keep = np.ones(nb, bool)
    for a, b in spans:
        keep[a:b] = False
    v = np.concatenate([tracks[t][keep & np.isfinite(tracks[t])] for t in tracks])
    if len(v):
        lo, hi = np.percentile(v, [0.5, 99.5])
        ax.set_ylim(lo - 0.1 * (hi - lo) - 1e-3, hi + 0.35 * (hi - lo) + 1e-3)
    sk = sct[keep & np.isfinite(sct)]
    if len(sk):
        a_, b_ = np.percentile(sk, [0.5, 99.5])
        ax2.set_ylim(a_ - 0.05 * (b_ - a_) - 1e-3, b_ + 0.05 * (b_ - a_) + 1e-3)
    ax.set_ylabel('NCP num fluctuation', fontsize=15)
    ax2.set_ylabel('Condensability', fontsize=15, rotation=-90, va='bottom')
    leg = ax.legend(loc='upper left', framealpha=1, fontsize=8, ncol=1 if len(tracks) <= 5 else 2)
    for h in leg.legend_handles:
        h.set_alpha(1.0); h.set_linewidth(2.0)
    ax.set_title(f'{series} {chrom}', fontsize=15)
    fig.tight_layout(); fig.savefig(out / 'ed2d.png', dpi=200); plt.close(fig)
    R['d'] = dict(chrom=chrom, mask=how, masked_bins=int((~keep).sum()), bins=int(nb),
                  score_point=int(sup), spearman_input_vs_score=_spearman(tracks[0], sct))


def panel_e(pos, seq, out, R, plt):
    pos = pos[(pos - 85 >= 0) & (pos + 86 <= len(seq))]
    if not len(pos):
        return False
    w = seq[pos[:, None] + np.arange(-85, 86)[None, :]]
    w = w[np.isin(w, np.frombuffer(b'ACGT', np.uint8)).all(axis=1)]
    comp = np.zeros(256, np.uint8)
    for a, b in zip(b'ACGT', b'TGCA'):
        comp[a] = b
    W, S = np.frombuffer(b'AT', np.uint8), np.frombuffer(b'CG', np.uint8)
    ww, ss, n = np.zeros(146), np.zeros(146), 0
    for i in (-1, 0, 1):                                       # dyad and +-1 bp, both strands (Park's notebook)
        sub = w[:, 85 + i - 73: 85 + i + 74]
        for q in (sub, comp[sub[:, ::-1]]):
            iw, is_ = np.isin(q, W), np.isin(q, S)
            ww += (iw[:, :-1] & iw[:, 1:]).sum(axis=0)
            ss += (is_[:, :-1] & is_[:, 1:]).sum(axis=0)
            n += len(q)
    if not n:
        return False
    ww, ss = ww / n, ss / n
    fig, ax = plt.subplots(figsize=(2.6, 2))
    ax.plot(ww, 'r', lw=2)
    ax2 = ax.twinx()
    ax2.plot(ss, 'b', lw=2)
    for k in range(-6, 7, 2):
        ax.axvline(73 + 10 * k, color='k', linestyle='--', alpha=0.25)
    ax.set_xticks([73 + 10 * k for k in range(-7, 8, 2)], [str(10 * k) for k in range(-7, 8, 2)], fontsize=6)
    ax.set_xlabel('Super Helical Location', fontsize=8)
    ax.set_ylabel('AA/AT/TA/TT frequency', fontsize=8, color='r')
    ax2.set_ylabel('CC/CG/GC/GG frequency', fontsize=8, color='b', rotation=-90, va='bottom')
    ax.tick_params('y', labelsize=5, colors='r'); ax2.tick_params('y', labelsize=5, colors='b')
    ax.set_title('Dinucleotide frequency', fontsize=8)
    fig.tight_layout(); fig.savefig(out / 'ed2e.png', dpi=200); plt.close(fig)
    R['e'] = dict(nucleosomes=int(len(w)), ww_period_bp=period(ww), ss_period_bp=period(ss),
                  ww=np.round(ww, 5).tolist(), ss=np.round(ss, 5).tolist())
    return True


def period(x, lo=5, hi=20):
    """Dominant period (bp) between lo and hi of a profile, after removing its 11 bp moving mean (the envelope
    across the nucleosome), from the zero-padded Fourier transform of the middle (edges of the mean dropped)."""
    x = np.asarray(x, float)
    r = (x - np.convolve(x, np.ones(11) / 11, 'same'))[5:-5]
    f = np.abs(np.fft.rfft(r - r.mean(), 4096))
    freq = np.fft.rfftfreq(4096)
    ok = (freq >= 1 / hi) & (freq <= 1 / lo)
    return float(1 / freq[ok][np.argmax(f[ok])])


def nuc_num(ncov_files, frac, mscale, chrom):
    """Per-nucleosome NCP number on chrom; total over all the files (NCPnum_ver2.py)."""
    tabs = [_read(p) for p in ncov_files]
    total = sum(t.iloc[:, 2].sum() for t in tabs)
    x = pd.concat(tabs)
    x = x[x.iloc[:, 0] == chrom].copy()
    x['num'] = np.round(MNUM * mscale * frac * x.iloc[:, 2] / max(total, 1))
    return x.iloc[:, [0, 1]].set_axis(['Chromosome', 'Position'], axis=1).assign(num=x['num'].to_numpy())


def panel_f(m, out, R, plt):
    X, Y, Z = m['num_in'].to_numpy(), m['num_sup'].to_numpy(), m['score'].to_numpy()
    fig, ax = plt.subplots(figsize=(6, 4))
    sca = ax.scatter(X + 1, Y + 1, c=Z, s=1, cmap='viridis', rasterized=True)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('NCP number (input)'); ax.set_ylabel('NCP number (sup)')
    cb = fig.colorbar(sca)
    cb.set_label('Condensability', rotation=-90, va='bottom')
    r_in, r_sup = _spearman(X, Z), _spearman(Y, Z)
    xm, ym = np.median(X + 1), np.median(Y + 1)
    y_lo, y_hi = np.percentile(Y + 1, [0.5, 99.9])
    x_lo, x_hi = np.percentile(X + 1, [1, 99.5])
    ax.annotate('', xy=(xm, y_hi), xytext=(xm, max(y_lo, 1.5)), arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
    ax.set_ylim(top=y_hi * 30)
    ax.text(xm, y_hi * 4, 'Corr = %.2f' % r_sup, color='red', fontweight='bold', ha='center')
    ax.annotate('', xy=(x_hi, ym), xytext=(x_lo, ym), arrowprops=dict(arrowstyle='->', color='magenta', lw=1.5))
    ax.text(x_hi * 1.15, ym, 'Corr = %.2f' % r_in, color='magenta', fontweight='bold', va='center')
    fig.tight_layout(); fig.savefig(out / 'ed2f.png', dpi=200); plt.close(fig)
    R['f'] = dict(spearman_score_input=r_in, spearman_score_sup=r_sup, nucleosomes=int(len(m)))


def _boxpair(ax, ax2, left, right, names, rot):
    pos = np.arange(len(names))
    for a, data, x, col in ((ax, left, pos - 0.2, 'tab:blue'), (ax2, right, pos + 0.2, 'tab:red')):
        bp = a.boxplot(data, notch=False, sym='', positions=x, widths=0.3, patch_artist=True, manage_ticks=False,
                       medianprops=dict(color='#ff7f0e'))
        for p in bp['boxes']:
            p.set_facecolor(col)
    ax.set_xticks(pos, names, rotation=rot, ha='right', rotation_mode='anchor')
    ax.set_xlim(-0.5, len(names) - 0.5)
    ax.tick_params('y', colors='tab:blue'); ax2.tick_params('y', colors='tab:red')


def panel_g(m, order, out, R, plt):
    import matplotlib.patches as mpatches
    inp = [m.loc[m.state == k, 'num_in'].to_numpy() for k, _ in order]
    sup = [m.loc[m.state == k, 'num_sup'].to_numpy() for k, _ in order]
    fig, ax = plt.subplots(figsize=(max(5, 0.33 * len(order) + 1), 3))
    ax2 = ax.twinx()
    _boxpair(ax, ax2, inp, sup, [v for _, v in order], 75)
    ax.set_ylabel('NCP number', color='tab:blue')
    ax2.set_ylabel('NCP number', color='tab:red', rotation=-90, va='bottom')
    ax.set_title('Estimated NCP number over chromatin states', fontsize=9)
    ax.legend(handles=[mpatches.Patch(facecolor='tab:blue', edgecolor='k', label='Input control'),
                       mpatches.Patch(facecolor='tab:red', edgecolor='k', label='Supernatant')],
              loc='upper right', fontsize=7)
    fig.tight_layout(); fig.savefig(out / 'ed2g.png', dpi=200); plt.close(fig)
    pair = lambda xs: [[cohen_abs(a, b) if i != j else np.nan for j, b in enumerate(xs)] for i, a in enumerate(xs)]
    dsup, din = np.array(pair(sup), float), np.array(pair(inp), float)
    R['g'] = dict(states=[v for _, v in order], n=[int(len(x)) for x in inp],
                  cohen_d_sup_max=float(np.nanmax(dsup)) if len(order) > 1 else None,
                  cohen_d_input_max=float(np.nanmax(din)) if len(order) > 1 else None)


def panel_h(m, win, order, out, R, plt):
    import matplotlib.patches as mpatches
    pk = [m.loc[m.state == k, 'score'].to_numpy() for k, _ in order]
    wn = [win.loc[win.state == k, 'score'].to_numpy() for k, _ in order]
    keep = [i for i in range(len(order)) if len(pk[i]) and len(wn[i])]
    pk, wn, order = [pk[i] for i in keep], [wn[i] for i in keep], [order[i] for i in keep]
    if not order:
        return False
    fig, ax = plt.subplots(figsize=(max(5, 0.33 * len(order) + 1), 3))
    ax2 = ax.twinx()
    _boxpair(ax, ax2, pk, wn, [v for _, v in order], 45)
    allp, allw = np.concatenate(pk), np.concatenate(wn)
    lo = min(np.percentile(allp, 1), np.percentile(allw, 1))
    hi = max(np.percentile(allp, 99.5), np.percentile(allw, 99.5))
    ax.set_ylim(lo - 0.3, hi + 1.4); ax2.set_ylim(lo - 0.3, hi + 1.4)
    ax.set_ylabel('Condensability(A.U.)', color='tab:blue')
    ax2.set_ylabel('Condensability(A.U.)', color='tab:red', rotation=-90, va='bottom')
    ax.set_title('NCP peak calling VS Sliding window', fontsize=9)
    stats = []
    for i, (a, b) in enumerate(zip(pk, wn)):
        d = cohen_abs(a, b)
        stats.append(dict(state=order[i][1], cohen_d=d, p=welch_p(a, b), n_peaks=int(len(a)), n_windows=int(len(b))))
        q = lambda x: np.percentile(x, 75) + 1.5 * (np.percentile(x, 75) - np.percentile(x, 25))
        y = min(max(q(a), q(b)), hi + 0.2)
        ax.plot([i - 0.2, i - 0.2, i + 0.2, i + 0.2], [y + 0.05, y + 0.15, y + 0.15, y + 0.05], color='red', lw=0.8)
        ax.text(i, y + 0.18, 'n.s.' if d < 0.1 else 'd=%.2f' % d, color='red', ha='center', fontsize=5)
    ax.text(0.98, 0.04, "n.s. Cohen's d < 0.1", transform=ax.transAxes, ha='right', fontsize=6)
    ax.legend(handles=[mpatches.Patch(facecolor='tab:blue', edgecolor='k', label='NCP peaks'),
                       mpatches.Patch(facecolor='tab:red', edgecolor='k', label='Sliding window')],
              loc='upper left', fontsize=7)
    fig.tight_layout(); fig.savefig(out / 'ed2h.png', dpi=200); plt.close(fig)
    R['h'] = dict(windows=int(len(win)), per_state=stats,
                  cohen_d_max=float(max(s['cohen_d'] for s in stats)))
    return True


def panel_i(m, order, out, R, plt):
    data = [m.loc[m.state == k, 'score'].to_numpy() for k, _ in order]
    n = len(data)
    up, lo = np.full((n, n), np.nan), np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i < j:
                up[i, j] = cohen_abs(data[i], data[j])
            elif i > j:
                lo[i, j] = -np.log10(welch_p(data[i], data[j]) + 1e-10)
    size = max(6.0, 0.42 * n + 2.0)
    fs = 7 if n <= 10 else 5.5 if n <= 15 else 4.5
    fig, ax = plt.subplots(figsize=(size + 1.6, size))
    im1 = ax.imshow(up, cmap='Reds')
    im2 = ax.imshow(lo, cmap='Blues', vmin=0, vmax=5)
    for i in range(n):
        for j in range(n):
            if i < j:
                ax.text(j, i, str(round(up[i, j], 1)), ha='center', va='center', color='tab:blue', fontsize=fs,
                        fontweight='bold')
            elif i > j:
                ax.text(j, i, str(round(lo[i, j], 1)), ha='center', va='center', color='tab:red', fontsize=fs - 1,
                        fontweight='bold')
    names = [v for _, v in order]
    ax.set_xticks(range(n), names, rotation=45, ha='right', rotation_mode='anchor')
    ax.set_yticks(range(n), names)
    c1 = fig.colorbar(im1, ax=ax, shrink=0.6, pad=0.02); c1.set_label("Cohen's d", rotation=-90, va='bottom')
    c2 = fig.colorbar(im2, ax=ax, shrink=0.6, pad=0.04); c2.set_label('-log10 p-value', rotation=-90, va='bottom')
    fig.tight_layout(); fig.savefig(out / 'ed2i.png', dpi=200); plt.close(fig)
    R['i'] = dict(states=names, cohen_d=np.round(up, 3).tolist(), neglog10p=np.round(lo, 2).tolist())


def window_scores(win_input, win_point, frac):
    """Condensability of the sliding windows: the per-nucleosome score (+1 pseudocount, windows with input, totals
    over them, minus ln f)."""
    c = _read(win_input).set_axis(['Chromosome', 'Start', 'End', 'v'], axis=1).set_index(['Start', 'End'])['v']
    t = _read(win_point).set_axis(['Chromosome', 'Start', 'End', 'v'], axis=1).set_index(['Start', 'End'])['v']
    t = t.reindex(c.index, fill_value=0)
    ok = c > 0
    c, t = c[ok], t[ok]
    sc = -np.log(((t + 1) / (t + 1).sum()) / ((c + 1) / (c + 1).sum())) - np.log(frac)
    return pd.DataFrame({'Start': sc.index.get_level_values(0), 'score': sc.to_numpy()})


# ---------------------------------------------------------------- driver

def run(series, samples, bins, bin_size, titration_path, fasta, chrom, outdir, agent='', mscale=1.0,
        ncov_input=(), ncov_point=(), nuc_scores=None, states=None, states_key=None, win_input=None,
        win_point=None, cytoband=None):
    """Draw the panels this series' data allow into outdir; return the summary (also written as ed2.json)."""
    mpl, plt = _plt()
    import pysam
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    tit = titration.read(titration_path)
    samples = [json.loads(Path(p).read_text()) for p in samples]
    R = {'series': series, 'chrom': chrom, 'panels': []}
    with pysam.FastaFile(str(fasta)) as fa:
        seq = np.frombuffer(fa.fetch(chrom).upper().encode(), dtype=np.uint8)
    panel_c(samples, tit, agent, out, R, plt, mpl); R['panels'].append('c')
    b10 = bins_10kb(bins, bin_size)
    sup = max(b10)
    bands = annotations.read_bands(cytoband, chrom) if cytoband else []
    panel_d(series, b10, tit, agent, mscale, chrom, seq, bands, sup, out, R, plt, mpl); R['panels'].append('d')
    if nuc_scores and ncov_input and ncov_point:
        sc = _read(nuc_scores)
        col = next(c for c in sc.columns[2:] if tnum_of(c) == sup)
        sc = sc[['Chromosome', 'Position', col]].rename(columns={col: 'score'}).dropna()
        if panel_e(sc['Position'].to_numpy(np.int64), seq, out, R, plt):
            R['panels'].append('e')
        a = nuc_num(ncov_input, tit.frac[0], mscale, chrom)
        b = nuc_num(ncov_point, tit.frac[sup], mscale, chrom)
        m = a.merge(b, on=['Chromosome', 'Position'], suffixes=('_in', '_sup')).merge(sc, on=['Chromosome', 'Position'])
        m = m.sort_values('Position')
        if len(m) > 2:
            panel_f(m, out, R, plt); R['panels'].append('f')
        if states and len(m) > 2:
            st = annotations.read_states(states, chrom)
            m = m.assign(state=categorize(m['Position'].to_numpy(), st))
            order = [(k, v) for k, v in annotations.state_order(states_key, {n for _, _, n in st})
                     if (m.state == k).sum() >= 2]
            R['states'] = dict(file=str(states), built_in=states_key, used=[v for _, v in order])
            if order:
                panel_g(m, order, out, R, plt); R['panels'].append('g')
                if win_input and win_point:
                    w = window_scores(win_input, win_point, tit.frac[sup])
                    w['state'] = categorize(w['Start'].to_numpy(), st)     # by window start, as Park's notebook
                    if panel_h(m, w, order, out, R, plt):
                        R['panels'].append('h')
                if len(order) > 1:
                    panel_i(m, order, out, R, plt); R['panels'].append('i')
    (out / 'ed2.json').write_text(json.dumps(R, indent=1, default=float))
    return R

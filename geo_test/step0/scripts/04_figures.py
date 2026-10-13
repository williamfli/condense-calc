"""Step 0 figures, from the outputs of 01 and 03 (and the manifest). Env: condrepro; login node is fine.

  figures/01_inputs_alignment.png     read pairs per run; bowtie2 outcome per sample
  figures/02_filters_fragments.png    Park's read filters; fragment lengths; MAPQ of kept fragments
  figures/03_bin_counts.png           fragments per 10 kb bin; chr1 counts and low-MAPQ share
  figures/04_scores_vs_geo.png        chr1 score track (ours, GEO, GEO rep1); scatter; differences
  figures/05_differences_by_mapq.png  |difference| by low-MAPQ class and by depth; seed test
  figures/06_genome_map.png           where the large differences are, with the low-MAPQ share
  figures/07_scale.png                ours-GEO and seed differences next to replicate differences
  figures/08_bowtie2_versions.png     each alignment of the same reads (bowtie2 version, seed) vs GEO
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from config import RESULTS, FIGURES, BIG, MANIFEST, CHROMS, LOWQ_LABELS

FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
C = {'ours': '#1f77b4', 'geo': '#d62728', 'rep1': '#7f7f7f', 'seed': '#2ca02c',
     'input': '#4c72b0', 'sup': '#dd8452'}
J = pd.read_csv(BIG / '03_bins_joined.tsv.gz', sep='\t')
J['lowq_class'] = pd.Categorical(J.lowq_class, LOWQ_LABELS, ordered=True)
J['Chromosome'] = pd.Categorical(J.Chromosome, CHROMS, ordered=True)
J = J.sort_values(['Chromosome', 'Start']).reset_index(drop=True)
mb = 1e6


def save(f, name):
    f.tight_layout()
    f.savefig(FIGURES / name, dpi=150)
    plt.close(f)
    print('wrote', name)


# ------------------------------------------------------------------ 01 inputs and alignment
man = pd.read_csv(MANIFEST, sep='\t')
man = man[man.step == 0].copy()
man['tnum'] = man.sample_title.str.split(',').str[3].astype(int)
bt = pd.read_csv(RESULTS / '03_bowtie2.tsv', sep='\t')
f, ax = plt.subplots(1, 2, figsize=(10, 3.4))
lab = [f'{r.run}\n{"input" if r.tnum == 0 else "ID 6"}' for r in man.itertuples()]
ax[0].bar(lab, man.read_pairs / 1e6, color=[C['input'] if t == 0 else C['sup'] for t in man.tnum])
for i, v in enumerate(man.read_pairs):
    ax[0].text(i, v / 1e6, f'{v / 1e6:.2f} M', ha='center', va='bottom')
ax[0].set_ylabel('read pairs (millions)')
ax[0].set_title('a  ENA runs (2 sequencing runs per library)')
b0 = bt[bt.alignment == 'ours'].set_index('sample')
parts = [('concordant_unique', 'concordant, 1 hit', '#4c72b0'),
         ('concordant_multi', 'concordant, >1 hit', '#c44e52'),
         ('not_concordant', 'not concordant', '#bbbbbb')]
x = np.arange(len(b0))
bottom = np.zeros(len(b0))
for col, name, color in parts:
    v = b0[col] / b0.pairs * 100
    ax[1].bar(x, v, bottom=bottom, color=color, label=name, width=0.6)
    for xi, (bi, vi) in enumerate(zip(bottom, v)):
        if vi > 3:
            ax[1].text(xi, bi + vi / 2, f'{vi:.1f}%', ha='center', va='center', color='w')
    bottom += v
ax[1].set_xticks(x, ['input (ID 0)', 'supernatant (ID 6)'])
ax[1].set_ylabel('% of read pairs')
ax[1].legend(frameon=False, fontsize=8, ncol=3, loc='upper center', bbox_to_anchor=(.5, -.1))
ax[1].set_title('b  bowtie2 (--trim-to 50, hg38 no-alt)')
save(f, '01_inputs_alignment.png')

# ------------------------------------------------------------------ 02 filters, lengths, MAPQ
casc = pd.read_csv(RESULTS / '01_filter_cascade.tsv', sep='\t', index_col=0)
tl = pd.read_csv(RESULTS / '01_tlen_hist.tsv', sep='\t', index_col=0)
mq = pd.read_csv(RESULTS / '01_mapq_hist.tsv', sep='\t', index_col=0)
f, ax = plt.subplots(1, 3, figsize=(13, 3.6))
y = np.arange(len(casc))
for k, (col, cc) in enumerate(zip(casc.columns, (C['input'], C['sup']))):
    ax[0].barh(y + (k - .5) * .38, casc[col] / 1e6, height=.38, color=cc, label=col)
    for yi, v in zip(y, casc[col]):
        ax[0].text(v / 1e6, yi + (k - .5) * .38, f' {v / 1e6:.2f}', va='center', fontsize=7)
ax[0].set_yticks(y, casc.index)
ax[0].invert_yaxis()
ax[0].set_xlabel('reads left (millions)')
ax[0].legend(frameon=False, fontsize=8, loc='lower right')
ax[0].set_title("a  Park's read filters, in order")
for col, cc in zip(tl.columns, (C['input'], C['sup'])):
    ax[1].plot(tl.index[:-1], tl[col].iloc[:-1] / tl[col].sum() * 100, color=cc, label=col, lw=1)
ax[1].axvspan(126, 184, color='0.9', zorder=0)
ax[1].text(155, ax[1].get_ylim()[1] * .95, 'mono-\nnucleosome\n126-184', ha='center', va='top', fontsize=7)
ax[1].set_xlim(0, 600)
ax[1].set_xlabel('fragment length |TLEN| (bp)')
ax[1].set_ylabel('% of kept fragments')
ax[1].set_title('b  fragment lengths (all are counted)')
w = 0.42
for k, (col, cc) in enumerate(zip(mq.columns, (C['input'], C['sup']))):
    ax[2].bar(mq.index + (k - .5) * w, mq[col] / mq[col].sum() * 100, width=w, color=cc, label=col)
ax[2].set_yscale('log')
ax[2].axvline(9.5, color='k', lw=.8, ls='--')
low = mq.loc[:9].sum() / mq.sum() * 100
ax[2].text(10, ax[2].get_ylim()[1] * .5 if ax[2].get_ylim()[1] > 1 else 50,
           f' MAPQ < 10:\n input {low.iloc[0]:.1f}%\n ID 6 {low.iloc[1]:.1f}%', fontsize=8, va='top')
ax[2].set_xlabel('MAPQ of the forward read')
ax[2].set_ylabel('% of kept fragments')
ax[2].set_title('c  mapping quality (no filter applied)')
save(f, '02_filters_fragments.png')

# ------------------------------------------------------------------ 03 counts per bin
f = plt.figure(figsize=(13, 5.6))
gs = f.add_gridspec(2, 3)
a = f.add_subplot(gs[:, 0])
bins = np.arange(0, 151, 2)
for col, cc, nm in (('input', C['input'], 'input (ID 0)'), ('sup', C['sup'], 'supernatant (ID 6)')):
    v = J[col].dropna()
    a.hist(v.clip(upper=150), bins=bins, histtype='step', color=cc, lw=1.2,
           label=f'{nm}: median {v.median():.0f}')
a.set_xlabel('fragments per 10 kb bin (150 = 150 or more)')
a.set_ylabel('bins')
a.legend(frameon=False, fontsize=8)
a.set_title('a  depth per bin')
c1 = J[J.Chromosome == 'chr1']
b = f.add_subplot(gs[0, 1:])
b.plot(c1.Start / mb, c1.input, color=C['input'], lw=.5, label='input')
b.plot(c1.Start / mb, c1.sup, color=C['sup'], lw=.5, label='ID 6')
b.set_ylim(0, np.nanpercentile(c1.input, 99.5) * 1.3)
b.set_ylabel('fragments / 10 kb')
b.legend(frameon=False, fontsize=8, ncol=2, loc='upper right')
b.set_title('b  chr1: fragments per bin')
c = f.add_subplot(gs[1, 1:], sharex=b)
c.fill_between(c1.Start / mb, c1.lowmapq_frac * 100, color='#c44e52', lw=0, step='mid')
c.set_ylabel('% MAPQ < 10')
c.set_xlabel('chr1 (Mb)')
c.set_title('c  chr1: share of fragments with MAPQ < 10 (input + ID 6)')
save(f, '03_bin_counts.png')

# ------------------------------------------------------------------ 04 scores vs GEO
f = plt.figure(figsize=(13, 7.5))
gs = f.add_gridspec(3, 3, height_ratios=[1, 1, 1.25])
for k, (lo_mb, hi_mb, what) in enumerate(((30, 40, 'uniquely mapping'), (143, 153, 'next to the centromere'))):
    a = f.add_subplot(gs[0, :2] if k == 0 else gs[0, 2:])
    w = c1[(c1.Start >= lo_mb * mb) & (c1.Start < hi_mb * mb)]
    a.plot(w.Start / mb, w.geo_rep1, color=C['rep1'], lw=.6, label='GEO replicate 1 (other reads)')
    a.plot(w.Start / mb, w.geo, color=C['geo'], lw=1.4, label='GEO replicate 2')
    a.plot(w.Start / mb, w.ours, color=C['ours'], lw=.8, label='ours (replicate 2 reads)')
    a.set_xlim(lo_mb, hi_mb)
    a.set_ylabel('score (ID 6)')
    a.set_xlabel('chr1 (Mb)')
    a.set_title(f'{"a" if k == 0 else "a2"}  chr1:{lo_mb}-{hi_mb} Mb, {what}')
    if k == 0:
        a.set_ylim(-0.2, None)
        a.legend(frameon=False, fontsize=8, ncol=3, loc='lower left')
b = f.add_subplot(gs[1, :])
b.plot(c1.Start / mb, c1.d, color='k', lw=.5)
b.axhline(0, color='0.6', lw=.5)
b.set_ylabel('ours - GEO')
b.set_xlabel('chr1 (Mb)')
b.set_title('b  difference along chr1 (same reads)')
ok = J.ours.notna() & J.geo.notna()
lo, hi = np.nanpercentile(J.geo[ok], [.05, 99.95])
for k, (xcol, ycol, nm, xl, yl) in enumerate((
        ('geo', 'ours', 'c  ours vs GEO, same reads', 'GEO replicate 2', 'ours'),
        ('ours', 'seed1', 'd  bowtie2 --seed 1 vs --seed 0, same reads', 'ours (seed 0)', 'ours (seed 1)'),
        ('geo_rep1', 'geo', 'e  GEO replicate 2 vs replicate 1', 'GEO replicate 1', 'GEO replicate 2'))):
    ax = f.add_subplot(gs[2, k])
    m = J[xcol].notna() & J[ycol].notna()
    hb = ax.hexbin(J[xcol][m], J[ycol][m], gridsize=90, extent=(lo, hi, lo, hi), norm=LogNorm(), cmap='viridis',
                   mincnt=1, linewidths=0)
    ax.plot([lo, hi], [lo, hi], color='w', lw=.6, ls='--')
    r = np.corrcoef(J[xcol][m], J[ycol][m])[0, 1]
    ax.text(.03, .95, f'r = {r:.3f}\nn = {m.sum():,}', transform=ax.transAxes, va='top', color='k',
            bbox=dict(fc='w', ec='none', alpha=.7))
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(nm, fontsize=9)
f.colorbar(hb, ax=ax, label='bins')
save(f, '04_scores_vs_geo.png')

# ------------------------------------------------------------------ 05 differences by MAPQ / depth
bym = pd.read_csv(RESULTS / '03_by_lowmapq.tsv', sep='\t')
f, ax = plt.subplots(1, 3, figsize=(13, 3.8))
data = [J.loc[J.lowq_class == c, 'd'].dropna().abs() for c in LOWQ_LABELS]
seedd = [J.loc[J.lowq_class == c, 'd_seed'].dropna().abs() for c in LOWQ_LABELS]
pos = np.arange(len(LOWQ_LABELS))
for dd, off, cc, nm in ((data, -.2, C['geo'], 'ours vs GEO'), (seedd, .2, C['seed'], 'seed 1 vs seed 0')):
    bp = ax[0].boxplot(dd, positions=pos + off, widths=.35, showfliers=False, patch_artist=True,
                       whis=(5, 95), medianprops=dict(color='k'))
    for p in bp['boxes']:
        p.set(facecolor=cc, alpha=.6)
    ax[0].plot([], [], color=cc, lw=6, alpha=.6, label=nm)
ax[0].set_xticks(pos, LOWQ_LABELS)
ax[0].set_xlabel('share of the bin\'s fragments with MAPQ < 10')
ax[0].set_ylabel('|difference| (box 25-75%, whiskers 5-95%)')
ax[0].legend(frameon=False, fontsize=8)
ax[0].set_title('a  |difference| by multi-mapper share')
g = bym[bym.lowq_class.isin(LOWQ_LABELS)]
for comp, cc, mk in (('ours vs GEO', C['geo'], 'o'), ('seed 1 vs seed 0', C['seed'], 's')):
    s = g[g.comparison == comp]
    ax[1].plot(pos, s.pearson, marker=mk, color=cc, label=comp)
ax2 = ax[1].twinx()
ax2.bar(pos, g[g.comparison == 'ours vs GEO'].share_of_bins * 100, color='0.85', zorder=0)
ax2.set_ylabel('% of bins (grey bars)')
ax[1].set_zorder(ax2.get_zorder() + 1)
ax[1].patch.set_visible(False)
ax[1].set_xticks(pos, LOWQ_LABELS)
ax[1].set_ylabel('Pearson r within the class')
ax[1].set_xlabel('share with MAPQ < 10')
ax[1].legend(frameon=False, fontsize=8, loc='lower left')
ax[1].set_title('b  agreement by multi-mapper share')
byd = pd.read_csv(RESULTS / '03_by_depth.tsv', sep='\t')
ax[2].plot(np.arange(len(byd)), byd.median_abs_diff, 'o-', color=C['geo'], label='median')
ax[2].plot(np.arange(len(byd)), byd.p95_abs_diff, 's--', color=C['geo'], label='95th percentile')
ax[2].set_xticks(np.arange(len(byd)), byd.sup_decile.str.replace(r'\(([-\d.]+), ([\d.]+)\]', lambda m:
                 f'{max(0, float(m.group(1))):.0f}-{float(m.group(2)):.0f}', regex=True), rotation=45, fontsize=7)
ax[2].set_xlabel('ID 6 fragments in the bin (deciles)')
ax[2].set_ylabel('|ours - GEO|')
ax[2].legend(frameon=False, fontsize=8)
ax[2].set_title('c  |difference| by depth')
save(f, '05_differences_by_mapq.png')

# ------------------------------------------------------------------ 06 genome map
f, ax = plt.subplots(len(CHROMS), 1, figsize=(13, 11), sharex=True)
lens = J.groupby('Chromosome').End.max()
for a, ch in zip(ax, CHROMS):
    g = J[J.Chromosome == ch]
    a.fill_between(g.Start / mb, g.lowmapq_frac, color='#c44e52', alpha=.35, lw=0, step='mid')
    big = g[g.d.abs() > 0.5]
    a.scatter(big.Start / mb, np.full(len(big), 1.05), s=4, color='k', marker='|', lw=.8)
    gap = g[g.geo.isna() & g.ours.isna()]
    a.scatter(gap.Start / mb, np.full(len(gap), -0.1), s=2, color='0.6', marker='|', lw=.5)
    a.set_xlim(0, lens.max() / mb)
    a.set_ylim(-0.2, 1.15)
    a.set_yticks([])
    a.spines['left'].set_visible(False)
    a.text(-2, .45, ch, ha='right', va='center', fontsize=7)
    a.hlines(-.1, 0, lens[ch] / mb, color='0.85', lw=.5, zorder=0)
ax[-1].set_xlabel('position (Mb)')
ax[0].set_title('Red: share of fragments with MAPQ < 10 (0-1). Black ticks: bins with |ours - GEO| > 0.5. '
                'Grey ticks: bins in neither table (no input reads).', fontsize=9)
save(f, '06_genome_map.png')

# ------------------------------------------------------------------ 07 scale
f, ax = plt.subplots(1, 2, figsize=(11, 3.6))
edges = np.linspace(-1.5, 1.5, 241)
series = (('ours - GEO (same reads)', J.d, C['geo']),
          ('seed 1 - seed 0 (same reads)', J.d_seed, C['seed']),
          ('GEO rep 2 - rep 1 (other reads)', J.geo - J.geo_rep1, C['rep1']))
for nm, v, cc in series:
    v = v.dropna()
    ax[0].hist(v.clip(-1.5, 1.5), bins=edges, histtype='step', color=cc, lw=1.2, label=f'{nm}: sd {v.std():.3f}')
ax[0].set_yscale('log')
ax[0].set_xlabel('difference in score (clipped to ±1.5)')
ax[0].set_ylabel('bins')
ax[0].legend(frameon=False, fontsize=8, loc='upper left')
ax[0].set_title('a  size of each kind of difference')
uniq = J.lowmapq_frac == 0
for nm, v, cc in series:
    q = [v[m].dropna().abs().quantile(.95) for m in (uniq, ~uniq & J.lowmapq_frac.notna())]
    ax[1].plot([0, 1], q, 'o-', color=cc, label=nm)
ax[1].set_xticks([0, 1], ['no MAPQ < 10 fragments', 'some MAPQ < 10 fragments'])
ax[1].set_xlim(-.4, 1.4)
ax[1].set_ylabel('95th percentile of |difference|')
ax[1].legend(frameon=False, fontsize=8)
ax[1].set_title('b  uniquely mapped bins vs the rest')
save(f, '07_scale.png')

# ------------------------------------------------------------------ 08 alignments of the same reads
av = pd.read_csv(RESULTS / '03_alignments_vs_geo.tsv', sep='\t')
av['label'] = [f'{r.bowtie2}\nseed {r.seed}' + ('' if r.index == 'NCBI' else f'\n{r.index} index')
               + ('\n(pipeline)' if r.alignment == 'ours' else '') for r in av.itertuples()]
f, ax = plt.subplots(1, 3, figsize=(13, 3.8))
x = np.arange(len(av))
ax[0].bar(x, av.same_counts_vs_geo * 100, color='#4c72b0')
for xi, v in zip(x, av.same_counts_vs_geo * 100):
    ax[0].text(xi, v, f'{v:.1f}%', ha='center', va='bottom', fontsize=8)
ax[0].set_ylabel('% of bins with the same counts as GEO')
ax[0].set_title('a  bins identical to GEO')
w = .38
ax[1].bar(x - w / 2, av.r_vs_geo_unique, width=w, color='#4c72b0', label='bins with no MAPQ < 10')
ax[1].bar(x + w / 2, av.r_vs_geo_gt50pct, width=w, color='#c44e52', label='bins > 50% MAPQ < 10')
ax[1].set_ylabel('Pearson r with GEO')
ax[1].set_ylim(0.5, 1.08)
ax[1].legend(frameon=False, fontsize=8, ncol=2, loc='upper center')
ax[1].set_title('b  agreement with GEO by MAPQ class')
ax[2].bar(x, av.same_counts_vs_ours.fillna(1) * 100, color='#2ca02c')
for xi, v in zip(x, av.same_counts_vs_ours.fillna(1) * 100):
    ax[2].text(xi, v, f'{v:.1f}%', ha='center', va='bottom', fontsize=8)
ax[2].set_ylabel('% of bins with the same counts as the pipeline')
ax[2].set_title('c  bins identical to the pipeline run')
for a in ax:
    a.set_xticks(x, av.label, fontsize=8)
    a.set_xlabel('bowtie2 version and seed' + (' (NCBI index)' if (av.index == 'NCBI').all() else ''))
save(f, '08_bowtie2_versions.png')

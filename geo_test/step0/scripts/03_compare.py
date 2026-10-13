"""Step 0, table level: compare our 10 kb score table with GEO's, and find where and why they differ.

Inputs
  RUN/tables/H1_NCP_PEG_2rep_10kb_score.gtab.gz       ours (pipeline, bowtie2 seed 0)
  SEED1/tables/H1_NCP_PEG_2rep_10kb_score.gtab.gz     ours, realigned with --seed 1 (02)
  GEO/GSE252941_H1_NCP_PEG_2rep_10kb_score.gtab.txt.gz  Park's table for the same reads
  GEO/GSE252941_H1_NCP_PEG_1rep_10kb_score.gtab.txt.gz  Park's biological replicate 1 (scale)
  RUN/bin/*_10kb_bin.gtab.gz                           our fragment counts per bin
  RUN/bt2_<version>/tables/...                         ours, realigned with older bowtie2 (05a)
  RUN/idx_4dn/tables/...                               ours, realigned to 4DN's bowtie2 index (06)
  BIG/01_bins_10kb_mapq.tsv.gz                         low-MAPQ fragments per bin (01)
  RUN/logs/align/*.log, <realignment>/logs/*.log       bowtie2 summaries

Outputs (results/)
  03_bowtie2.tsv          alignment summary per sample and seed
  03_comparisons.tsv      pairwise comparisons of score tables
  03_alignments_vs_geo.tsv  each alignment (bowtie2 version / seed) against GEO, overall and by MAPQ
  03_offset.tsv           the constant shift between ours and GEO, and the share of bins that match
  03_by_lowmapq.tsv       ours vs GEO and seed 1 vs seed 0, by low-MAPQ class
  03_by_depth.tsv         ours vs GEO by supernatant fragments per bin (deciles)
  03_by_chrom.tsv         ours vs GEO per chromosome
  03_row_sets.tsv         bins present in only one table, and why
  03_top_discordant.tsv   the 40 bins with the largest |ours - GEO|
  BIG/03_bins_joined.tsv.gz  everything per bin, for 04_figures.py
Env: condrepro. Login node is fine (~1 min, ~3 GB).
"""
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from config import (RUN, SEED1, BT2_VERSIONS, IDX4DN, GEO, BIG, RESULTS, COL, KEY, CHROMS, LOWQ_EDGES, LOWQ_LABELS)

RESULTS.mkdir(exist_ok=True)


def table(path, col=COL, name=None):
    t = pd.read_csv(path, sep='\t', na_values='NA').set_index(KEY)
    return t[col].rename(name or col)


# ------------------------------------------------------------------ bowtie2 summaries
def bowtie2(path):
    txt = open(path).read()
    num = lambda pat: int(re.search(pat, txt).group(1))
    return {'pairs': num(r'(\d+) reads; of these'),
            'concordant_unique': num(r'(\d+) \([\d.]+%\) aligned concordantly exactly 1 time'),
            'concordant_multi': num(r'(\d+) \([\d.]+%\) aligned concordantly >1 times'),
            'not_concordant': num(r'(\d+) \([\d.]+%\) aligned concordantly 0 times'),
            'overall_rate_pct': float(re.search(r'([\d.]+)% overall alignment rate', txt).group(1))}


# alignments of the same reads: name -> (bowtie2 version, seed, index, log dir, table dir)
ALIGNS = {'ours': ('2.5.5', 0, 'NCBI', RUN / 'logs/align', RUN / 'tables'),
          'seed1': ('2.5.5', 1, 'NCBI', SEED1 / 'logs', SEED1 / 'tables')}
for v in BT2_VERSIONS:
    if (RUN / f'bt2_{v}/tables').exists():
        ALIGNS[f'bt2_{v}'] = (v, 0, 'NCBI', RUN / f'bt2_{v}/logs', RUN / f'bt2_{v}/tables')
if (IDX4DN / 'tables').exists():
    ALIGNS['idx_4dn'] = ('2.5.5', 0, '4DN', IDX4DN / 'logs', IDX4DN / 'tables')
rows = []
for name, (ver, seed, index, logs, _) in ALIGNS.items():
    for tnum in (0, 6):
        rows.append({'alignment': name, 'bowtie2': ver, 'seed': seed, 'index': index,
                     'sample': f'H1_NCP_PEG_{tnum}_2rep',
                     **bowtie2(logs / f'H1_NCP_PEG_{tnum}_2rep.log')})
bt = pd.DataFrame(rows)
bt['multi_share_of_concordant'] = bt.concordant_multi / (bt.concordant_unique + bt.concordant_multi)
bt.to_csv(RESULTS / '03_bowtie2.tsv', sep='\t', index=False, float_format='%.4f')
print(bt.to_string())

# ------------------------------------------------------------------ one table per bin
alns = [table(t / 'H1_NCP_PEG_2rep_10kb_score.gtab.gz', name=n) for n, (*_, t) in ALIGNS.items()]
geo = table(GEO / 'GSE252941_H1_NCP_PEG_2rep_10kb_score.gtab.txt.gz', name='geo')
geo1 = table(GEO / 'GSE252941_H1_NCP_PEG_1rep_10kb_score.gtab.txt.gz', 'H1_NCP_PEG_6_1rep', 'geo_rep1')
cnt = {}
for tnum, nm in ((0, 'input'), (6, 'sup')):
    cnt[nm] = table(RUN / f'bin/H1_NCP_PEG_{tnum}_2rep_10kb_bin.gtab.gz', f'H1_NCP_PEG_{tnum}_2rep', nm)
mq = pd.read_csv(BIG / '01_bins_10kb_mapq.tsv.gz', sep='\t')
# matched on (chromosome, start): the last bin's End is the chromosome end in 01, maybe not elsewhere
mq = mq.pivot_table(index=['Chromosome', 'Start'], columns='sample', values=['fragments', 'mapq_lt10'],
                    aggfunc='sum')
lowq = mq['mapq_lt10'].sum(axis=1) / mq['fragments'].sum(axis=1).clip(lower=1)

J = pd.concat(alns + [geo, geo1, cnt['input'], cnt['sup']], axis=1)
J = J[J.index.get_level_values(0).isin(CHROMS)]
J = J.iloc[np.lexsort((J.index.get_level_values('Start'), J.index.get_level_values(0).map(CHROMS.index)))]
cs = J.index.droplevel('End')
J['lowmapq_frac'] = lowq.reindex(cs).to_numpy()
check = (mq['fragments']['H1_NCP_PEG_0_2rep'].reindex(cs).fillna(0).to_numpy() - J.input.fillna(0)).abs().max()
print(f'01 fragment counts vs bincount input column: max |difference| = {check}')
J['lowq_class'] = pd.cut(J.lowmapq_frac, LOWQ_EDGES, labels=LOWQ_LABELS)
J['d'] = J.ours - J.geo
J['d_seed'] = J.seed1 - J.ours
J.reset_index().to_csv(BIG / '03_bins_joined.tsv.gz', sep='\t', index=False, float_format='%.6g')


def compare(a, b, label):
    x, y = J[a], J[b]
    both = x.notna() & y.notna()
    d = (x - y)[both]
    return {'comparison': label, 'bins_both': int(both.sum()),
            f'only_first': int((x.notna() & y.isna()).sum()), 'only_second': int((y.notna() & x.isna()).sum()),
            'pearson': np.corrcoef(x[both], y[both])[0, 1], 'spearman': spearmanr(x[both], y[both])[0],
            'mean_diff': d.mean(), 'sd_diff': d.std(), 'median_abs_diff': d.abs().median(),
            'p95_abs_diff': d.abs().quantile(0.95), 'frac_abs_lt_0.05': (d.abs() < 0.05).mean(),
            'frac_abs_gt_0.5': (d.abs() > 0.5).mean()}


cmp = pd.DataFrame([compare('ours', 'geo', 'ours vs GEO (same reads)'),
                    compare('seed1', 'ours', 'ours seed 1 vs ours seed 0 (same reads)'),
                    compare('ours', 'geo_rep1', 'ours vs GEO replicate 1 (other reads)'),
                    compare('geo', 'geo_rep1', 'GEO replicate 2 vs GEO replicate 1')]
                   + [compare(n, r, f'bowtie2 {ALIGNS[n][0]}, {ALIGNS[n][2]} index vs {lab} (same reads)')
                      for n in ALIGNS if n.startswith(('bt2_', 'idx_'))
                      for r, lab in (('geo', 'GEO'), ('ours', 'ours (bowtie2 2.5.5)'))])
cmp.to_csv(RESULTS / '03_comparisons.tsv', sep='\t', index=False, float_format='%.4f')
print(cmp.round(4).to_string())

# ------------------------------------------------------------------ the constant shift
# If a bin has the same counts in both tables, d = ln(S_ours/S_geo) - ln(C_ours/C_geo): the same
# constant for every such bin (S, C: totals of s+1 and c+1). Scores are printed to 5 decimals.
d = J.d.dropna()
h = np.round(d / 1e-5).astype(int).value_counts()
delta = h.index[0] * 1e-5
same = (d - delta).abs() <= 2e-5
inp = J.input[J.input > 0]
S_ours, C_ours = (J.sup[J.input > 0] + 1).sum(), (inp + 1).sum()
off = pd.DataFrame([{'mode_of_diff': delta, 'bins_at_mode_pm_2e-5': int(same.sum()),
                     'share_of_bins_at_mode': same.mean(), 'mean_diff': d.mean(),
                     'S_ours_sum_sup_plus1': S_ours, 'C_ours_sum_input_plus1': C_ours}])
off.to_csv(RESULTS / '03_offset.tsv', sep='\t', index=False, float_format='%.6g')
J['same_counts'] = ((J.d - delta).abs() <= 2e-5)


def at_mode(x):
    """Share of bins whose difference is the table-wide constant (same counts, see above)."""
    x = x.dropna()
    m = np.round(x / 1e-5).astype(int).value_counts().index[0] * 1e-5
    return ((x - m).abs() <= 2e-5).mean(), m


av = []
for name, (ver, seed, index, *_) in ALIGNS.items():
    row = {'alignment': name, 'bowtie2': ver, 'seed': seed, 'index': index}
    for ref in ('geo', 'ours'):
        if ref == name:
            continue
        dd = J[name] - J[ref]
        ok = dd.notna()
        row[f'r_vs_{ref}'] = np.corrcoef(J[name][ok], J[ref][ok])[0, 1]
        row[f'same_counts_vs_{ref}'], row[f'shift_vs_{ref}'] = at_mode(dd)
        row[f'frac_abs_gt_0.5_vs_{ref}'] = (dd[ok].abs() > .5).mean()
        for lab_c, m in (('unique', J.lowmapq_frac == 0), ('gt50pct', J.lowmapq_frac > .5)):
            k = ok & m
            row[f'r_vs_{ref}_{lab_c}'] = np.corrcoef(J[name][k], J[ref][k])[0, 1]
    av.append(row)
av = pd.DataFrame(av)
av.to_csv(RESULTS / '03_alignments_vs_geo.tsv', sep='\t', index=False, float_format='%.5g')
print(av.T.to_string())
print(off.to_string())
J.reset_index().to_csv(BIG / '03_bins_joined.tsv.gz', sep='\t', index=False, float_format='%.6g')


# ------------------------------------------------------------------ stratified
def strat(g, dcol):
    x = g[dcol].dropna()
    a, b = (g.ours, g.geo) if dcol == 'd' else (g.seed1, g.ours)
    ok = a.notna() & b.notna()
    return pd.Series({'bins': len(x), 'share_of_bins': np.nan,
                      'pearson': np.corrcoef(a[ok], b[ok])[0, 1] if ok.sum() > 2 else np.nan,
                      'median_abs_diff': x.abs().median(), 'p95_abs_diff': x.abs().quantile(.95),
                      'frac_abs_lt_0.05': (x.abs() < .05).mean(), 'frac_abs_gt_0.5': (x.abs() > .5).mean(),
                      'same_counts_share': g.loc[x.index, 'same_counts'].mean() if dcol == 'd' else np.nan})


out = []
for dcol, lab in (('d', 'ours vs GEO'), ('d_seed', 'seed 1 vs seed 0')):
    s = J.groupby('lowq_class', observed=True).apply(lambda g: strat(g, dcol))
    s['share_of_bins'] = s.bins / s.bins.sum()
    out.append(s.assign(comparison=lab).reset_index())
    cum = []
    for lab_c, mx in (('<= 2%', .02), ('<= 10%', .1), ('all', 1.0)):
        g = J[J.lowmapq_frac <= mx]
        cum.append(strat(g, dcol).rename(f'cumulative {lab_c}'))
    c = pd.DataFrame(cum)
    c['share_of_bins'] = c.bins / J[dcol].notna().sum()
    out.append(c.assign(comparison=lab).rename_axis('lowq_class').reset_index())
bym = pd.concat(out)[['comparison', 'lowq_class', 'bins', 'share_of_bins', 'pearson', 'median_abs_diff',
                      'p95_abs_diff', 'frac_abs_lt_0.05', 'frac_abs_gt_0.5', 'same_counts_share']]
bym.to_csv(RESULTS / '03_by_lowmapq.tsv', sep='\t', index=False, float_format='%.4f')
print(bym.round(3).to_string())

J['sup_decile'] = pd.qcut(J.sup, 10, duplicates='drop')
byd = J.groupby('sup_decile', observed=True).apply(lambda g: strat(g, 'd')).drop(columns='share_of_bins')
byd.to_csv(RESULTS / '03_by_depth.tsv', sep='\t', float_format='%.4f')
byc = J.groupby(level=0).apply(lambda g: strat(g, 'd')).drop(columns='share_of_bins')
byc['median_lowmapq_frac'] = J.groupby(level=0).lowmapq_frac.median()
byc = byc.loc[CHROMS]
byc.to_csv(RESULTS / '03_by_chrom.tsv', sep='\t', float_format='%.4f')
print(byc.round(3).to_string())

# ------------------------------------------------------------------ rows in one table only
rs = []
for name, m in (('only GEO', J.geo.notna() & J.ours.isna()), ('only ours', J.ours.notna() & J.geo.isna())):
    g = J[m]
    rs.append({'set': name, 'bins': int(m.sum()), 'our_input_0': int((g.input == 0).sum()),
               'our_input_1': int((g.input == 1).sum()), 'our_input_ge2': int((g.input >= 2).sum()),
               'on_chrY': int((g.index.get_level_values(0) == 'chrY').sum()),
               'median_lowmapq_frac': g.lowmapq_frac.median()})
pd.DataFrame(rs).to_csv(RESULTS / '03_row_sets.tsv', sep='\t', index=False, float_format='%.3f')
print(pd.DataFrame(rs).to_string())

top = J.dropna(subset=['d']).reindex(J.d.abs().dropna().sort_values(ascending=False).index[:40])
top[['ours', 'geo', 'd', 'seed1', 'd_seed', 'input', 'sup', 'lowmapq_frac']].to_csv(
    RESULTS / '03_top_discordant.tsv', sep='\t', float_format='%.4f')

# ------------------------------------------------------------------ uniquely mapped bins that still differ
# Size of the difference (after removing the constant shift) against the change one fragment makes.
u = J[(J.lowmapq_frac == 0) & J.d.notna()]
v = u[~u.same_counts]
x = (v.d - delta).abs()
one = np.log((v.sup.clip(lower=1) + 1) / v.sup.clip(lower=1))      # one more/less supernatant fragment
res = pd.DataFrame([{'unique_bins': len(u), 'differing': len(v), 'share_differing': len(v) / len(u),
                     'median_abs_diff_minus_shift': x.median(), 'p95_abs_diff_minus_shift': x.quantile(.95),
                     'median_one_fragment_change': one.median(),
                     'share_within_1.5x_one_fragment': (x <= 1.5 * one).mean()}])
res.to_csv(RESULTS / '03_unique_residual.tsv', sep='\t', index=False, float_format='%.4f')
print(res.T.to_string())

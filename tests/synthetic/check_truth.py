"""Do the pipeline's tables recover the simulated truth? (tests/simulate.py)

    python tests/synthetic/check_truth.py OUT       (OUT = the folder given to tests/synthetic/run.sh)

Truth per nucleosome i at concentration c: score = -ln p_i(c), p_i = 1/(1+exp(rate (c - chalf_i))).
Binned truth: -ln of the bin's mean p (scores measure the bin's surviving fraction).
Prints Pearson r between measured and true values for 10 kb / 1 kb scores, per-nucleosome
scores (dyads matched within 20 bp) and c-half.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(sys.argv[1])
truth = pd.read_csv(OUT / 'input/truth/nucleosomes.tsv', sep='\t')
tit = pd.read_csv(OUT / 'input/titration/SIM_NCP_sp.tsv', sep='\t').set_index('tnum')
R = OUT / 'run/tables'


def p_of(conc, chalf, rate):
    return 1 / (1 + np.exp(rate * (conc - chalf)))


def binned(kb, path, label_fmt):
    t = pd.read_csv(path, sep='\t')
    res = []
    for col in t.columns[3:]:
        tnum = int(col.split('_')[3])
        conc = tit.loc[tnum, 'conc']
        tr = truth.assign(p=p_of(conc, truth.chalf, truth.rate), b=truth.dyad // (kb * 1000))
        mean_p = tr.groupby(['chrom', 'b'])['p'].mean()
        key = list(zip(t.Chromosome, t.Start // (kb * 1000)))
        true = -np.log(mean_p.reindex(key).to_numpy())
        ok = np.isfinite(true) & np.isfinite(t[col].to_numpy(dtype=float))
        res.append((tnum, np.corrcoef(t[col].to_numpy(dtype=float)[ok], true[ok])[0, 1]))
    print(f'{label_fmt}: ' + ', '.join(f'tnum {n}: r={r:.3f}' for n, r in res))


binned(10, R / 'SIM_NCP_sp_1rep_10kb_score.gtab.gz', 'scores, 10 kb bins (shallow)')
binned(1, R / 'SIM_NCP_sp_1rep_deep_1kb_score.gtab.gz', 'scores, 1 kb bins (deep)  ')

# per-nucleosome scores: match called dyads to true dyads
rows = []
for f in sorted((R / 'SIM_NCP_sp_1rep_deep_score').glob('*.gtab.gz')):
    rows.append(pd.read_csv(f, sep='\t'))
ncp = pd.concat(rows)
called = 0
for col in ncp.columns[2:]:
    tnum = int(col.split('_')[3])
    conc = tit.loc[tnum, 'conc']
    rs = []
    for chrom, g in ncp.groupby('Chromosome'):
        tr = truth[truth.chrom == chrom].sort_values('dyad')
        idx = np.searchsorted(tr.dyad.to_numpy(), g.Position.to_numpy())
        idx = np.clip(idx, 1, len(tr) - 1)
        left, right = tr.dyad.to_numpy()[idx - 1], tr.dyad.to_numpy()[idx]
        near = np.where(np.abs(g.Position - left) <= np.abs(g.Position - right), idx - 1, idx)
        d = np.abs(tr.dyad.to_numpy()[near] - g.Position.to_numpy())
        m = (d <= 20) & np.isfinite(g[col].to_numpy(dtype=float))       # NA: a dyad without input coverage
        true = -np.log(p_of(conc, tr.chalf.to_numpy()[near], tr.rate.to_numpy()[near]))
        rs.append(pd.DataFrame({'meas': g[col].to_numpy(dtype=float)[m], 'true': true[m]}))
        called = max(called, m.sum())
    d = pd.concat(rs)
    print(f'scores, per nucleosome, tnum {tnum}: r={np.corrcoef(d.meas, d.true)[0, 1]:.3f} '
          f'({len(d)} of {len(ncp)} called dyads within 20 bp of a true dyad; {len(truth)} true nucleosomes)')

# c-half per 10 kb bin vs the bin's mean true c-half
fit = pd.read_csv(R / 'SIM_NCP_sp_1rep_10kb_4PL.gtab.gz', sep='\t', na_values='NA')
tb = truth.assign(b=truth.dyad // 10000).groupby(['chrom', 'b'])['chalf'].mean()
true = tb.reindex(list(zip(fit.Chromosome, fit.Start // 10000))).to_numpy()
ok = fit['C-half'].notna().to_numpy() & np.isfinite(true)
print(f'c-half, 10 kb bins: r={np.corrcoef(fit["C-half"][ok], true[ok])[0, 1]:.3f} '
      f'(fitted {ok.sum()} of {len(fit)} bins; mean |error| {np.mean(np.abs(fit["C-half"][ok] - true[ok])):.3f} mM)')

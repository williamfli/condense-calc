"""Compare a pipeline table with the matching GEO GSE252941 table.

    python geo_test/compare_geo.py OURS.gtab.gz GEO.gtab.txt.gz

Rows are matched on (chromosome, start, end) or (chromosome, position). For each column present
in both: rows only in ours / only in GEO, Pearson and Spearman r, mean difference (ours - GEO),
median |difference|, and the share of rows within 0.05.
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def load(path):
    t = pd.read_csv(path, sep='\t', na_values='NA')
    key = list(t.columns[:2]) if t.columns[1] == 'Position' else list(t.columns[:3])
    return t.set_index(key), key


ours, key = load(sys.argv[1])
geo, _ = load(sys.argv[2])
print(f'ours: {len(ours):,} rows   GEO: {len(geo):,} rows   key: {", ".join(key)}')
print(f'{"column":28s} {"both":>9s} {"only ours":>9s} {"only GEO":>9s} {"pearson":>8s} {"spearman":>8s}'
      f' {"mean diff":>9s} {"med |d|":>8s} {"|d|<0.05":>8s}')
for col in [c for c in ours.columns if c in geo.columns]:
    a, b = ours[col].dropna(), geo[col].dropna()
    both = a.index.intersection(b.index)
    x, y = a.loc[both].to_numpy(float), b.loc[both].to_numpy(float)
    d = x - y
    print(f'{col:28s} {len(both):9,d} {len(a.index.difference(b.index)):9,d} {len(b.index.difference(a.index)):9,d}'
          f' {np.corrcoef(x, y)[0, 1]:8.4f} {spearmanr(x, y)[0]:8.4f} {d.mean():9.4f} {np.median(np.abs(d)):8.4f}'
          f' {np.mean(np.abs(d) < 0.05):8.1%}')
missing = [c for c in geo.columns if c not in ours.columns]
if missing:
    print('columns only in GEO:', ', '.join(missing))

"""Compare our outputs with the original scripts' outputs on the synthetic data, stage by stage.

    python tests/equivalence/compare.py OUT      (OUT = .../condense-seq-pipeline/synthetic)

For every file pair: same header? same rows (keys) in the same order? how many value cells are
identical as text, and the largest numeric difference. Writes OUT/equiv/comparison.tsv.
"""
import math
import sys
from pathlib import Path

from condseq import gtab

OUT = Path(sys.argv[1])
O, R, E = OUT / 'equiv/original', OUT / 'run', OUT / 'equiv/ours'
CHR = ['chr1', 'chr2', 'chrX']
INP = 'SIM_NCP_sp_0_1rep_deep'
names = sorted(p.name[:-4] for p in (R / 'bam').glob('*.bam'))
deep = [n for n in names if n.endswith('_deep')]

pairs = []
for n in names:
    for kb in ('1kb', '10kb'):
        pairs.append(('bincount', f'{n}_{kb}_bin.gtab.gz', R / 'bin'))
        if '_sp_0_' not in n:
            pairs.append(('score (bins)', f'{n}_{kb}_score.gtab.gz', R / 'score'))
for n in deep:
    for c in CHR:
        pairs.append(('coverage', f'{n}_{c}_cov.gtab.gz', E))
        pairs.append(('ncov', f'{n}_{c}_Ncov.gtab.gz', R / 'areas'))
        pairs.append(('bindata', f'{n}_{c}_171win25step_Bdata.gtab.gz', R / 'areas'))
        if '_sp_0_' not in n:
            pairs.append(('score (nucleosomes)', f'{n}_{c}_score.gtab.gz', R / 'score'))
for c in CHR:
    pairs.append(('peaks', f'{INP}_{c}_peak.gtab.gz', R / 'peaks'))
for t in range(6):
    pairs.append(('num', f'SIM_NCP_sp_{t}_1rep_10kb_num.gtab.gz', R / 'num'))
pairs.append(('fit (4PL)', 'SIM_NCP_sp_1rep_10kb_4PL.gtab.gz', R / 'tables'))


def num(x):
    try:
        return float(x)
    except ValueError:
        return x


def compare(a, b):
    ta, tb = gtab.read(a), gtab.read(b)
    nk = ta.nkey
    keys_a = [tuple(r[:nk]) for r in ta.rows]
    keys_b = [tuple(r[:nk]) for r in tb.rows]
    same_keys = keys_a == keys_b
    cells = same = 0
    maxdiff = 0.0
    if same_keys:
        for ra, rb in zip(ta.rows, tb.rows):
            for x, y in zip(ra[nk:], rb[nk:]):
                cells += 1
                same += x == y
                vx, vy = num(x), num(y)
                if isinstance(vx, float) and isinstance(vy, float):
                    if not (math.isnan(vx) and math.isnan(vy)):
                        maxdiff = max(maxdiff, abs(vx - vy))
                elif vx != vy:
                    maxdiff = math.inf
    return dict(header=ta.header == tb.header, rows_orig=len(ta.rows), rows_ours=len(tb.rows),
                same_rows=same_keys, cells=cells, identical_text=same, max_abs_diff=maxdiff)


lines = ['stage\tfile\theader\trows_orig\trows_ours\tsame_rows\tcells\tidentical_text\tmax_abs_diff']
summary = {}
for stage, fname, ours_dir in pairs:
    a, b = O / fname, ours_dir / fname
    if not a.exists() or not b.exists():
        lines.append(f'{stage}\t{fname}\tMISSING orig={a.exists()} ours={b.exists()}')
        summary.setdefault(stage, []).append(None)
        continue
    r = compare(a, b)
    lines.append('\t'.join([stage, fname] + [str(v) for v in r.values()]))
    summary.setdefault(stage, []).append(r)
(OUT / 'equiv/comparison.tsv').write_text('\n'.join(lines) + '\n')

print(f'{"stage":22s} files  rows-identical  text-identical cells        max |diff|')
for stage, rs in summary.items():
    ok = [r for r in rs if r]
    rows_ok = sum(r['same_rows'] and r['header'] for r in ok)
    cells = sum(r['cells'] for r in ok)
    same = sum(r['identical_text'] for r in ok)
    md = max((r['max_abs_diff'] for r in ok), default=float('nan'))
    miss = len(rs) - len(ok)
    print(f'{stage:22s} {len(rs):5d}  {rows_ok:5d}/{len(ok):<5d}    {same:>10d}/{cells:<10d}  {md:.3g}'
          + (f'   ({miss} missing)' if miss else ''))

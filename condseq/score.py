"""Scoring stages: counts -> condensability scores, molecule numbers, and column concatenation.

Ports (external/condense-seq/prepro_scripts, a8d2dd8):
  score()        NCPscore_ver4.py --numc   condensability = -ln(normalized sup / normalized input) - ln(f)
  numbers()      NCPnum_ver2.py            estimated nucleosome numbers per bin / dyad
  concatenate()  concat_gtab.py            outer join of per-titration files into one table
"""
import math

from . import gtab

INPUT_MOLECULES = 1.6e12            # NCPnum_ver2.py: input_mnum = mscale * 1.6e12
AGENT_MSCALE = {'sp': 1, 'spd': 1, 'CoH': 1, 'PEG': 1, 'Ca': 5, 'HP1a': 1, 'HP1bSUV': 1}  # batch_submit_ver7.py


def _keyed(table):
    """{row key: row} with keys that sort like the original's merge (chrom order, then positions)."""
    nk = table.nkey
    return {gtab.point_key(r[0], *[int(v) for v in r[1:nk]]): r for r in table.rows}


def score(pairs, survival_fraction=None, chroms=None, missing_as_nan=False):
    """Condensability scores of one titration point.

    pairs: list of (test_path, control_path) gtabs (one pair for a bin table; one per
           chromosome for per-nucleosome Ncov tables). Normalization totals run over all pairs.
    survival_fraction: soluble fraction f at this titration point (the --numc correction,
           used for GEO); None = no correction.
    Follows NCPscore_ver4.py exactly:
      - only rows present in the control are scored; a test row missing there is dropped;
        a control row missing in the test counts as 0
      - a row whose control values are all 0 is dropped (missing_as_nan=True writes it as NA instead,
        which is not what GEO has)
      - per column: +1 pseudocount on both sides; totals T, C over the rows with control > 0
      - score = -ln[((t+1)/T) / ((c+1)/C)] - ln(f), rounded to 5 decimals; control 0 -> NA
    Returns [(header, rows)] per pair, ready for gtab.write.
    """
    tables = []
    for tpath, cpath in pairs:
        t, c = gtab.read(tpath, chroms), gtab.read(cpath, chroms)
        assert t.header[:t.nkey] == c.header[:c.nkey], (tpath, cpath)
        tables.append((t, c))
    ncol = len(tables[0][0].data_columns())

    def merged(t, c):
        tk, ck = _keyed(t), _keyed(c)
        tcols = ccols = t.data_columns()     # the original reads the control at the test's column indices
        for key in sorted(ck):
            crow = ck[key]
            trow = tk.get(key)
            tvals = [float(trow[i]) for i in tcols] if trow is not None else [0.0] * ncol
            cvals = [float(crow[i]) for i in ccols]
            yield crow[:c.nkey], tvals, cvals

    # pass 1: totals
    T, C = [0.0] * ncol, [0.0] * ncol
    for t, c in tables:
        for _, tv, cv in merged(t, c):
            if sum(cv) <= 0:
                continue
            for i in range(ncol):
                if cv[i] > 0:
                    T[i] += tv[i] + 1
                    C[i] += cv[i] + 1
    offset = -math.log(float(survival_fraction) / float(1)) if survival_fraction is not None else 0

    # pass 2: scores
    out = []
    for t, c in tables:
        header = t.header[:t.nkey] + [t.header[i] for i in t.data_columns()]
        rows = []
        for key, tv, cv in merged(t, c):
            if sum(cv) <= 0:
                if missing_as_nan:
                    rows.append(list(key) + ['NA'] * ncol)
                continue
            s = []
            for i in range(ncol):
                if cv[i] <= 0:
                    s.append('NA')
                    continue
                ncount = float(tv[i] + 1) / T[i]
                ncontrol = float(cv[i] + 1) / C[i]
                val = -math.log(ncount / ncontrol) + offset
                s.append(gtab.py2_str(gtab.py2_round(val, 5)))
            rows.append(list(key) + s)
        out.append((header, rows))
    return out


def numbers(paths, survival_fraction, mscale=1, chroms=None):
    """Estimated nucleosome numbers of one titration point (NCPnum_ver2.py):
    num = round(mscale * 1.6e12 * f * value / total), total over all rows of all `paths`.
    Returns [(header, rows)] per path."""
    tables = [gtab.read(p, chroms) for p in paths]
    ncol = len(tables[0].data_columns())
    total = [0.0] * ncol
    for t in tables:
        for r in t.rows:
            for j, i in enumerate(t.data_columns()):
                total[j] += float(r[i])
    mnum = mscale * INPUT_MOLECULES * survival_fraction
    out = []
    for t in tables:
        cols = t.data_columns()
        header = t.header[:t.nkey] + [t.header[i] for i in cols]
        rows = [r[:t.nkey] + [int(gtab.py2_round(mnum * (float(r[i]) / total[j])))
                              for j, i in enumerate(cols)] for r in t.rows]
        out.append((header, rows))
    return out


def concatenate(paths, colnums=None):
    """Outer join of gtabs on their keys (concat_gtab.py). colnums: per file, the column
    indices to take (default: all data columns). Missing rows give 'NA'."""
    tables = [gtab.read(p) for p in paths]
    nk = tables[0].nkey
    if colnums is None:
        colnums = [t.data_columns() for t in tables]
    keyed = [_keyed(t) for t in tables]
    header = tables[0].header[:nk] + [t.header[i] for t, cols in zip(tables, colnums) for i in cols]
    rows = []
    for key in sorted(set().union(*keyed)):
        row = None
        vals = []
        for k, cols in zip(keyed, colnums):
            r = k.get(key)
            if r is None:
                vals += ['NA'] * len(cols)
            else:
                row = row or r[:nk]
                vals += [r[i] for i in cols]
        rows.append(row + vals)
    return header, rows

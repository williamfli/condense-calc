"""c-half: logistic fit of the surviving fraction against agent concentration, per bin.

Port of logistic_fit_ver3.py (already Python 3 in the original). Defaults are the GEO run's
(submit_scripts/logistic-submit, batch_submit_ver7.py): sigmoid model, differential evolution
(seed 3), top fixed at 1, bottom fixed at 0, rate in [0, 100], c-half in [min conc, max conc],
fits with R^2 < 0.5 reported as NA. Input: a concatenated table of nucleosome numbers
(score.numbers + score.concatenate) whose columns are the titration points, input included.
"""
import warnings

import numpy as np
from scipy.optimize import curve_fit, differential_evolution

from . import gtab
from .settings import DEFAULTS

FIT = DEFAULTS['fit']


def sigmoid(x, top, rate, chalf, bottom):
    return bottom + (top - bottom) / (1 + np.exp(rate * (x - chalf)))


def hill(x, top, rate, chalf, bottom):
    return bottom + (top - bottom) / (1.0 + (x / float(chalf)) ** rate)


MODELS = {'sigmoid': sigmoid, 'hill': hill}


def _sse(params, func, X, Y):
    warnings.filterwarnings('ignore')
    return np.sum((func(X, *params) - Y) ** 2)


def fit_table(path, titration, tnums, model=FIT['model'], method=FIT['method'], min_rsq=FIT['min_rsq'],
              top=(1, 1), bottom=(0, 0), rate=(0, 100), chalf=(None, None), chroms=None):
    """Fit every row of a num table. tnums: the titration ID of each data column, in order.
    Bounds given as None are set as in the original (from the data). Returns (header, rows)."""
    t = gtab.read(path, chroms)
    cols = t.data_columns()
    assert len(cols) == len(tnums), (len(cols), tnums)
    concs = [titration.conc[n] for n in tnums]
    order = [i for _, i in sorted((c, i) for i, c in enumerate(concs))]
    inp = order[0]
    assert concs[inp] == 0, 'the input column (concentration 0) is missing'
    X = np.asarray([concs[i] for i in order])
    func = MODELS[model]
    header = t.header[:t.nkey] + ['Model', 'Method', 'Top', 'Rate', 'C-half', 'Bottom', 'R-squared']
    rows = []
    for r in t.rows:
        if any(r[i] == 'NA' for i in cols):
            continue   # the original would crash on NA; concatenated num tables have none in practice
        nums = [int(r[i]) for i in cols]
        if nums[inp] <= 0:
            continue
        Y = np.asarray([nums[i] / nums[inp] for i in order], dtype=float)
        lo_top, hi_top = top[0] if top[0] is not None else max(Y) * 0.99, top[1] if top[1] is not None else max(Y) * 1.01
        lo_bot, hi_bot = bottom[0] if bottom[0] is not None else min(Y) * 0.99, bottom[1] if bottom[1] is not None else min(Y) * 1.01
        lo_ch, hi_ch = chalf[0] if chalf[0] is not None else min(X), chalf[1] if chalf[1] is not None else max(X)
        lo_r, hi_r = rate[0] if rate[0] is not None else 0.0, rate[1] if rate[1] is not None else 100.0
        bounds = [(lo_top, hi_top), (lo_r, hi_r), (lo_ch, hi_ch), (lo_bot, hi_bot)]
        if method == 'curve_fit':
            p0 = [max(Y), 1.0, float(np.mean(X)), min(Y)]
            try:
                p, _ = curve_fit(func, X, Y, p0, bounds=tuple(zip(*bounds)), method='dogbox')
                ok = True
            except Exception:
                ok = False
        else:
            res = differential_evolution(_sse, args=(func, X, Y), bounds=bounds, seed=3)
            p, ok = res.x, res.success
        if ok:
            resid = Y - func(X, *p)
            rsq = 1 - np.sum(resid ** 2) / np.sum((Y - np.mean(Y)) ** 2)
            ok = rsq >= min_rsq
        vals = [str(v) for v in (*p, rsq)] if ok else ['NA'] * 5
        rows.append(r[:t.nkey] + [model, method] + vals)
    return header, rows

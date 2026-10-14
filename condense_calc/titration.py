"""Titration tables: spermine (or other agent) concentration and soluble fraction per titration point.

The scores need, for every titration point, the fraction of nucleosomes left in the
supernatant (f, "survival fraction", measured by UV-Vis), and the c-half fit needs the
concentration. Park's own titration files were not published (GEO has none), so the
pipeline reads a simple table we write ourselves:

  tnum    conc    survival_fraction
  0       0       1.0
  1       0.158   0.93
  ...

tnum is the titration ID used in file and column names (0 = input). Park's 10-column
files (`<cell>_<sample>_<agent>_titration.csv`) are also accepted: concentration in
column 1, total survival fraction in column 8 (= third from last), tnum in the last column,
as NCPscore_ver4.py and logistic_fit_ver3.py read them.
"""
from dataclasses import dataclass


@dataclass
class Titration:
    conc: dict       # tnum -> concentration
    frac: dict       # tnum -> survival fraction


def read(path):
    with open(path) as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip()]
    head = lines[0].split()
    conc, frac = {}, {}
    if head[:3] == ['tnum', 'conc', 'survival_fraction']:
        for line in lines[1:]:
            t, c, f = line.split()[:3]
            conc[int(t)], frac[int(t)] = float(c), float(f)
    else:  # Park format, skip rows whose last column isn't an int (as the original does)
        for line in lines:
            cols = line.split('\t') if '\t' in line else line.split()
            try:
                t = int(cols[-1])
            except ValueError:
                continue
            conc[t], frac[t] = float(cols[0]), float(cols[-3])
    if 0 not in frac:
        raise ValueError(f'{path}: no row for tnum 0 (the input)')
    return Titration(conc, frac)

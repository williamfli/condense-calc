"""Counting stages: fragments -> per-bin counts, coverage, nucleosome peaks, coverage areas.

Ports (external/condense-seq/prepro_scripts, a8d2dd8):
  bincount()   bincount_ver3.py     reads per genomic bin (+ GC content, mean fragment length)
  coverage()   coverage_ver3.py     fragment coverage per bp (kept in memory, not written by default)
  peaks()      NCPpeak_ver3.py      nucleosome dyads, greedy non-overlapping selection
  ncpcov()     NCPcov_ver3.py       coverage area in a 171-bp window around each dyad
  bindata()    bindata_ver3.py      coverage area in sliding windows (171 bp every 25 bp)
Defaults are the values the GEO run used (submit_scripts/batch_submit_ver7.py and the *-submit
job scripts); each function says which.
"""
import functools
import math

import numpy as np
from sortedcontainers import SortedList

# Length windows used for GEO GSE252941 (bincount/coverage: none; NCPpeak: its defaults).
PEAK_MIN_LEN, PEAK_MAX_LEN = 126, 184
PEAK_NCP_LEN, PEAK_OVERLAP = 147, 40
NCP_WINDOW = 171                       # batch_submit_ver7.py: Nlen
SLIDING = (171, 25)                    # batch_submit_ver7.py: sliding_windows


# ---------------------------------------------------------------- bin counts

def bincount(frags, chrom_sizes, chroms, win_size, gc=None, tlen=True):
    """Reads per bin of `win_size` bp, for every bin of every chromosome in `chroms`.

    As bincount_ver3.py with `-w win_size --gc --tlen` and no length filter (the GEO run):
    one count per fragment at its centre (dyad_bincount), bins [i*w, (i+1)*w).
    Returns rows [chrom, start, end, count, (GC), (mean |tlen|)].
    `gc` is {chrom: array of GC fractions per bin} from gc_content(), or None.
    """
    rows = []
    for chrom in chroms:
        n = int(math.ceil(chrom_sizes[chrom] / win_size))
        counts = np.zeros(n, dtype=np.int64)
        tl = np.zeros(n, dtype=np.int64)
        f = frags.get(chrom)
        if f is not None and len(f):
            idx = f.dyad_bincount() // win_size
            # the original indexes a Python list: negative indices wrap, too-large ones crash
            if (idx >= n).any():
                raise IndexError(f'{chrom}: fragment centre beyond the last bin')
            idx = np.where(idx < 0, idx + n, idx)
            np.add.at(counts, idx, 1)
            np.add.at(tl, idx, np.abs(f.tlen))
        for i in range(n):
            row = [chrom, i * win_size, (i + 1) * win_size, int(counts[i])]
            if gc is not None:
                row.append('%f' % gc[chrom][i])
            if tlen:
                row.append('%f' % (tl[i] / counts[i] if counts[i] > 0 else 0))
            rows.append(row)
    return rows


def gc_content(fasta, chroms, win_size):
    """GC fraction of each bin, as bincount_ver3.py --gc computes it (upper-cased G+C over
    bin length, Ns counted in the length; the last partial bin uses its own length)."""
    import pysam
    out = {}
    with pysam.FastaFile(str(fasta)) as fa:
        for chrom in chroms:
            seq = np.frombuffer(fa.fetch(chrom).upper().encode(), dtype=np.uint8)
            is_gc = (seq == ord('G')) | (seq == ord('C'))
            n = int(math.ceil(len(seq) / win_size))
            sums = np.add.reduceat(is_gc.astype(np.int64), np.arange(0, len(seq), win_size)) if len(seq) else np.zeros(0)
            lens = np.minimum(win_size, len(seq) - np.arange(n) * win_size)
            out[chrom] = sums / lens
    return out


# ---------------------------------------------------------------- coverage

def coverage(frags, chrom_size):
    """Fragment coverage per bp of one chromosome (coverage_ver3.py; the GEO run used no
    length filter). Returns an int64 array long enough to hold every fragment end."""
    start, stop = frags.interval()
    n = max(chrom_size, int(stop.max()) + 1 if len(stop) else 0)
    if len(start) and start.min() < 0:
        raise ValueError('fragment starting before position 0')
    diff = np.zeros(n + 1, dtype=np.int64)
    np.add.at(diff, start, 1)
    np.add.at(diff, stop, -1)
    return np.cumsum(diff[:-1])


def coverage_rows(chrom, cov):
    """The rows coverage_ver3.py --skip writes: positions with non-zero coverage."""
    pos = np.flatnonzero(cov)
    return [[chrom, int(p), int(cov[p])] for p in pos]


# ---------------------------------------------------------------- peaks

def _tuple_cmp(a, b):
    # NCPpeak_ver3.py's comparator: never returns 0, so tied scores are ordered by the sort
    # algorithm's internals. Python 3 (3.6 and 3.12 tested) reproduces Python 2.7's order
    # exactly through cmp_to_key, so the same peaks are picked.
    return -1 if a[0] <= b[0] else 1


def peaks(frags, min_len=PEAK_MIN_LEN, max_len=PEAK_MAX_LEN, ncp_len=PEAK_NCP_LEN,
          overlap=PEAK_OVERLAP, fast=False):
    """Nucleosome dyads of one chromosome, as NCPpeak_ver3.py selects them.

    Dyad scores: mono-nucleosome fragments (126-184 bp by default), weight 1 at the centre
    (0.5/0.5 for even lengths). Positions are visited from the highest score down and kept
    unless a kept dyad lies closer than ncp_len - overlap (107 bp). No smoothing.
    fast=True breaks ties by position instead of reproducing Python 2's order (faster, but
    can pick different dyads among equal scores).
    Returns {position: score}.
    """
    f = frags.select(min_len, max_len)
    if not len(f):
        return {}
    pos, w = f.dyad_split()
    upos, inv = np.unique(pos, return_inverse=True)
    score = np.bincount(inv, weights=w)
    # Python 2 kept a position's score an int until a 0.5 was added; the text output shows it
    is_float = np.bincount(inv, weights=(w != 1)) > 0
    value = {int(p): (float(s) if fl else int(s)) for p, s, fl in zip(upos, score, is_float)}
    if fast:
        order = np.lexsort((upos, -score))
        cand = [(float(score[i]), int(upos[i])) for i in order]
    else:
        temp = [[float(s), int(p)] for s, p in zip(score, upos)]   # position order, as the original
        cand = sorted(temp, key=functools.cmp_to_key(_tuple_cmp), reverse=True)
    min_gap = ncp_len - overlap
    kept = SortedList()
    out = {}
    for s, p in cand:
        i = kept.bisect_left(p)
        if i > 0 and abs(kept[i - 1] - p) < min_gap:
            continue
        if i < len(kept) and abs(kept[i] - p) < min_gap:
            continue
        kept.add(p)
        out[p] = value[p]
    return out


# ---------------------------------------------------------------- coverage areas

def ncpcov(cov, dyads, window=NCP_WINDOW):
    """Coverage summed over [dyad - window//2, dyad + window//2] for each dyad
    (NCPcov_ver3.py --Nlen 171: 171 positions). Returns an int64 array aligned with dyads."""
    half = window // 2
    cs = np.concatenate([[0], np.cumsum(cov)])
    d = np.asarray(dyads, dtype=np.int64)
    lo = np.clip(d - half, 0, len(cov))
    hi = np.clip(d + half + 1, 0, len(cov))
    return cs[hi] - cs[lo]


def bindata(cov, chrom_size, size=SLIDING[0], step=SLIDING[1]):
    """Coverage summed over windows [k*step, k*step + size), clipped at the chromosome end
    (bindata_ver3.py --Bsize 171 --Bstep 25, sum mode). Like the original, only windows that
    contain covered positions are returned. Returns (starts, ends, sums)."""
    c = cov[:chrom_size]
    cs = np.concatenate([[0], np.cumsum(c)])
    k = np.arange(chrom_size // step + 1)
    st = k * step
    ed = np.minimum(st + size, chrom_size)
    sums = cs[np.minimum(ed, len(c))] - cs[np.minimum(st, len(c))]
    keep = sums > 0
    return st[keep], ed[keep], sums[keep]

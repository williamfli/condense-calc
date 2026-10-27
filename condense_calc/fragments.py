"""Read a BAM into per-chromosome fragment arrays, with the read filters of the original scripts.

Ports the read loop shared by bincount_ver3.py, coverage_ver3.py and NCPpeak_ver3.py
(`samtools view -F 0x10` + the checks after it). Each fragment is represented by its
forward-strand read, so a properly paired fragment is seen once.

Kept, as in the original:
  - forward-strand reads only (flag 0x10 unset)
  - mapped (0x4 unset), mate mapped (0x8 unset), proper pair (0x2 set), POS >= 1
  - edit distance NM <= mm_cutoff (a read without an NM tag is kept)
Not applied by the original, available here as opt-in options (default off = faithful):
  - min_mapq: MAPQ filter (the original's is commented out, so multi-mappers are kept)
  - dedup: drop fragments with the same (chromosome, start, template length)
Fragment-length windows differ by stage, so they are applied later (see select()).
"""
from array import array
from dataclasses import dataclass

import numpy as np
import pysam

from .settings import DEFAULTS

READS = DEFAULTS['reads']


@dataclass
class Fragments:
    pos: np.ndarray    # 0-based leftmost position of the forward read
    end: np.ndarray    # 0-based end of the forward read on the reference (pos + M + D)
    tlen: np.ndarray   # SAM TLEN of the forward read (negative if it is the right-hand read)

    def __len__(self):
        return len(self.pos)

    def select(self, min_len=0, max_len=None):
        """Fragments with min_len <= |tlen| <= max_len (the original's --min/--max)."""
        a = np.abs(self.tlen)
        keep = a >= min_len
        if max_len is not None:
            keep &= a <= max_len
        return Fragments(self.pos[keep], self.end[keep], self.tlen[keep])

    def interval(self):
        """Fragment [start, end) on the reference, as coverage_ver3.py computes it."""
        right = self.tlen <= 0
        start = np.where(right, self.end + self.tlen, self.pos)
        stop = np.where(right, self.end, self.pos + self.tlen)
        return start, stop

    def dyad_bincount(self):
        """Fragment centre as bincount_ver3.py assigns it (one position, weight 1).
        Even lengths go to the left of the two centre bases. Python 2 `/` on ints floors,
        which is `//` here, also for negative tlen."""
        t = self.tlen
        left = t > 0
        d_left = np.where(t % 2 != 0, self.pos + t // 2, self.pos + t // 2 - 1)
        d_right = self.end + t // 2 - 1
        return np.where(left, d_left, d_right)

    def dyad_split(self):
        """Fragment centre as NCPpeak_ver3.py assigns it: odd lengths give one position
        with weight 1; even lengths split 0.5/0.5 over the two centre bases.
        Returns (positions, weights), both 1-D, possibly longer than the input."""
        t = self.tlen
        left = t > 0
        odd = t % 2 != 0
        base = np.where(left, self.pos + t // 2, self.end + t // 2)
        # odd: left -> pos + t//2 ; right -> end + t//2 - 1
        p_odd = np.where(left, base, base - 1)[odd]
        # even: left -> (pos + t//2 - 1, pos + t//2) ; right -> (end + t//2 - 1, end + t//2)
        b_even = base[~odd]
        pos = np.concatenate([p_odd, b_even - 1, b_even])
        w = np.concatenate([np.ones(len(p_odd)), np.full(2 * len(b_even), 0.5)])
        return pos, w


def _ref_span(cigartuples):
    # the original adds M and D lengths only (bowtie2 emits no =/X/N)
    return sum(n for op, n in cigartuples if op in (0, 2))


def read_bam(path, chroms, mm_cutoff=READS['mm'], min_mapq=READS['min_mapq'], dedup=READS['dedup']):
    """Return {chrom: Fragments} for the chromosomes in `chroms` (None = all)."""
    want = set(chroms) if chroms else None
    acc = {}
    with pysam.AlignmentFile(str(path), 'rb') as bam:
        if want is not None and bam.has_index():
            # same reads as a full pass: each read is fetched once, from its own chromosome
            reads = (r for c in chroms if c in bam.references for r in bam.fetch(c))
        else:
            reads = bam.fetch(until_eof=True)
        for r in reads:
            flag = r.flag
            if flag & 0x10 or flag & 0x4 or flag & 0x8 or not flag & 0x2:
                continue
            ref = r.reference_name
            if ref is None or (want is not None and ref not in want):
                continue
            pos = r.reference_start
            if pos < 0:
                continue
            if r.has_tag('NM') and r.get_tag('NM') > mm_cutoff:
                continue
            if min_mapq and r.mapping_quality < min_mapq:
                continue
            cig = r.cigartuples or []
            if ref not in acc:
                acc[ref] = array('q')          # flat (pos, end, tlen) triples: 24 bytes per fragment
            acc[ref].extend((pos, pos + _ref_span(cig), r.template_length))
    out = {}
    for ref, buf in acc.items():
        a = np.frombuffer(buf, dtype=np.int64).reshape(-1, 3)
        if dedup and len(a):
            _, first = np.unique(a[:, [0, 2]], axis=0, return_index=True)
            a = a[np.sort(first)]
        out[ref] = Fragments(a[:, 0], a[:, 1], a[:, 2])
    return out

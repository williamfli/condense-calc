"""The gtab ("genomic table") format used by the Park 2025 scripts and GEO GSE252941.

Tab-separated, gzipped, with a header row. Two layouts:
  point:   Chromosome  Position  <value columns...>
  binned:  Chromosome  Start  End  <value columns...>   (0-based, half-open)
Some binned files end with extra GCcontent / Meantlen columns that are not data.
"""
import gzip
import math
from decimal import Decimal, ROUND_HALF_UP


def chrom_key(name):
    """Sort order of the original chr_cmp: numbered chromosomes first (numerically),
    then the rest by name (chrM, chrX, chrY). Python 2 put every int before every str."""
    suffix = name[3:] if name.startswith('chr') else name
    try:
        return (0, int(suffix), '')
    except ValueError:
        return (1, 0, suffix)


def sort_chroms(names):
    return sorted(names, key=chrom_key)


def point_key(chrom, *pos):
    """Row order used when the original scripts merge two sorted gtabs."""
    return (chrom_key(chrom),) + tuple(pos)


def py2_str(x):
    """str() of a number as Python 2 printed it: ints as-is, floats with %.12g plus '.0'
    when that looks like an int. Used so our text output matches the original byte for byte."""
    if isinstance(x, str):
        return x
    if isinstance(x, (int,)) and not isinstance(x, bool):
        return str(x)
    x = float(x)
    if math.isnan(x):
        return 'nan'
    if math.isinf(x):
        return 'inf' if x > 0 else '-inf'
    s = '%.12g' % x
    if not any(c in s for c in '.e'):
        s += '.0'
    return s


def py2_round(x, ndigits=0):
    """Python 2.7 round(): rounds the exact binary value, halves away from zero
    (Python 3 rounds halves to even). Decimal(x) is exact, so this matches 2.7 bit for bit."""
    q = Decimal(1).scaleb(-ndigits)
    return float(Decimal(x).quantize(q, rounding=ROUND_HALF_UP))


def open_text(path, mode='rt'):
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


class Table:
    """A gtab held in memory as rows of (key strings, value strings)."""

    def __init__(self, header, rows):
        self.header = header          # list of column names
        self.rows = rows              # list of lists of strings

    @property
    def layout(self):
        if self.header[1] == 'Position':
            return 'point'
        assert self.header[1] == 'Start' and self.header[2] == 'End', self.header[:3]
        return 'binned'

    @property
    def nkey(self):
        return 2 if self.layout == 'point' else 3

    def data_columns(self):
        """Indices of value columns: everything after the key, up to GCcontent if present."""
        end = self.header.index('GCcontent') if 'GCcontent' in self.header else len(self.header)
        return list(range(self.nkey, end))


def read(path, chroms=None):
    header, rows = None, []
    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            if header is None:
                header = cols
                continue
            if chroms and cols[0] not in chroms:
                continue
            rows.append(cols)
    return Table(header, rows)


def write(path, header, rows):
    with open_text(path, 'wt') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(py2_str(v) for v in row) + '\n')

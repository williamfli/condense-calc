"""Small hand-checkable tests of each stage. Run: pytest tests/test_units.py (in the cloudpoint env)."""
import math

import numpy as np

from cloudpoint import counts, gtab, score
from cloudpoint.fragments import Fragments


def F(rows):
    a = np.asarray(rows, dtype=np.int64).reshape(-1, 3)
    return Fragments(a[:, 0], a[:, 1], a[:, 2])


def test_dyad_rules():
    # left read, odd 147: centre at pos + 73; even 148: left of the two centre bases
    # right read (tlen < 0), odd -147: end + (-74) - 1
    f = F([[100, 150, 147], [100, 150, 148], [250, 300, -147]])
    assert list(f.dyad_bincount()) == [173, 173, 225]
    pos, w = f.dyad_split()
    assert sorted(zip(pos.tolist(), w.tolist())) == [(173, 0.5), (173, 1.0), (174, 0.5), (225, 1.0)]


def test_interval():
    f = F([[100, 150, 147], [250, 300, -147]])
    s, e = f.interval()
    assert list(s) == [100, 153] and list(e) == [247, 300]


def test_bincount():
    f = F([[100, 150, 147], [1000, 1050, 147], [2100, 2150, 147]])
    rows = counts.bincount({'chr1': f}, {'chr1': 2500}, ['chr1'], 1000, gc=None, tlen=True)
    assert [r[3] for r in rows] == [1, 1, 1]            # dyads 173, 1073, 2173
    assert [r[1:3] for r in rows] == [[0, 1000], [1000, 2000], [2000, 3000]]
    assert rows[0][4] == '147.000000'


def test_peaks_greedy():
    # dyad 1000 (3 reads) beats 1050 (2 reads, 50 bp away); 1200 (1 read) is 200 bp away, kept
    rows = [[1000 - 73, 1000, 147]] * 3 + [[1050 - 73, 1050, 147]] * 2 + [[1200 - 73, 1200, 147]]
    pk = counts.peaks(F(rows))
    assert pk == {1000: 3.0, 1200: 1.0}
    assert counts.peaks(F(rows), fast=True) == pk


def test_ncpcov_and_bindata():
    cov = np.zeros(1000, dtype=np.int64)
    cov[100:300] = 2
    # window [200-85, 200+85] = 171 positions, all inside the covered stretch
    assert counts.ncpcov(cov, [200])[0] == 2 * 171
    st, ed, sm = counts.bindata(cov, 1000, size=171, step=25)
    assert st[0] == 0 and sm[0] == 2 * 71                # [0,171) overlaps [100,171)
    assert (sm > 0).all()


def test_score_formula(tmp_path):
    t, c = tmp_path / 't.gtab.gz', tmp_path / 'c.gtab.gz'
    gtab.write(t, ['Chromosome', 'Start', 'End', 'T'], [['chr1', 0, 10, 5], ['chr1', 10, 20, 0], ['chr1', 30, 40, 9]])
    gtab.write(c, ['Chromosome', 'Start', 'End', 'C'], [['chr1', 0, 10, 10], ['chr1', 10, 20, 4], ['chr1', 20, 30, 0]])
    (header, rows), = score.score([(t, c)], survival_fraction=0.5)
    # rows scored: control rows with control > 0 (chr1:0, chr1:10); the test-only row is dropped
    T, C = (5 + 1) + (0 + 1), (10 + 1) + (4 + 1)
    want0 = -math.log(((5 + 1) / T) / ((10 + 1) / C)) - math.log(0.5)
    assert [r[1] for r in rows] == ['0', '10']
    assert float(rows[0][3]) == round(want0, 5)


def test_py2_formatting():
    assert gtab.py2_str(3) == '3' and gtab.py2_str(3.0) == '3.0' and gtab.py2_str(2.04069) == '2.04069'
    assert gtab.py2_round(0.125, 2) == 0.13          # Python 3 round() gives 0.12
    assert gtab.sort_chroms(['chrX', 'chr10', 'chr2', 'chrM', 'chr1']) == ['chr1', 'chr2', 'chr10', 'chrM', 'chrX']


def test_concatenate_stream_matches_memory(tmp_path):
    # sorted inputs take the streaming merge; it must give the in-memory outer join row for row
    hdr = ['Chromosome', 'Start', 'End']
    a = [['chr1', '0', '10', '1'], ['chr2', '0', '10', '2'], ['chr10', '0', '10', '3'], ['chrX', '0', '10', '4']]
    b = [['chr2', '0', '10', '5', '6', '0.4'], ['chr10', '0', '10', '7', '8', '0.5'], ['chrM', '0', '10', '9', '9', '0.6']]
    gtab.write(tmp_path / 'a.gtab.gz', hdr + ['a'], a)
    gtab.write(tmp_path / 'b.gtab.gz', hdr + ['b1', 'b2', 'GCcontent'], b)
    paths = [tmp_path / 'a.gtab.gz', tmp_path / 'b.gtab.gz']
    assert all(score._sorted_unique(p) for p in paths)
    h1, r1 = score._concatenate_stream(paths)
    h2, r2 = score._concatenate_memory(paths)
    r1, r2 = list(r1), list(r2)
    assert h1 == h2 == hdr + ['a', 'b1', 'b2'] and r1 == r2
    assert r1[1] == ['chr2', '0', '10', '2', '5', '6'] and r1[3] == ['chrM', '0', '10', 'NA', '9', '9']
    # an unsorted input is detected, so concatenate falls back to the in-memory join
    gtab.write(tmp_path / 'c.gtab.gz', hdr + ['c'], [a[2], a[0]])
    assert not score._sorted_unique(tmp_path / 'c.gtab.gz')
    h, r = score.concatenate([tmp_path / 'c.gtab.gz', paths[0]])
    assert [x[0] for x in r] == ['chr1', 'chr2', 'chr10', 'chrX']

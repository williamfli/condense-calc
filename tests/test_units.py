"""Small hand-checkable tests of each stage. Run: pytest tests/test_units.py (in the condense-calc env)."""
import math

import numpy as np
import pytest

from condense_calc import counts, gtab, score, sheet
from condense_calc.fragments import Fragments


@pytest.fixture(autouse=True)
def no_installs(tmp_path, monkeypatch):
    """Nothing is downloaded or installed: the default bowtie2 (2.4.1) would be, into ~/.cache. Tests that check
    an install mock subprocess.run themselves."""
    import subprocess
    from condense_calc import reference
    monkeypatch.setenv('CONDENSE_REFS', str(tmp_path / 'refs'))
    real = subprocess.run

    def run(cmd, *a, **k):
        if 'create' in cmd or any('get_reference' in str(c) for c in cmd):
            pytest.fail(f'a test would install: {cmd}')
        return real(cmd, *a, **k)
    monkeypatch.setattr(reference.subprocess, 'run', run)


def F(rows):
    a = np.asarray(rows, dtype=np.int64).reshape(-1, 3)
    return Fragments(a[:, 0], a[:, 1], a[:, 2])


def test_dyad_rules():
    # left read, odd 147: centre at pos + 73; even 148: left of the two centre bases
    # right read (tlen < 0), odd -147, fragment [153, 300): centre 226; Park: end + (-74) - 1 = 225
    # right read, even -148, fragment [152, 300): left centre base 225, split 225/226, in both modes
    f = F([[100, 150, 147], [100, 150, 148], [250, 300, -147], [252, 300, -148]])
    assert list(f.dyad_bincount()) == [173, 173, 226, 225]
    assert list(f.dyad_bincount(park=True)) == [173, 173, 225, 225]
    pos, w = f.dyad_split()
    assert sorted(zip(pos.tolist(), w.tolist())) == \
        [(173, 0.5), (173, 1.0), (174, 0.5), (225, 0.5), (226, 0.5), (226, 1.0)]
    pos, w = f.dyad_split(park=True)
    assert sorted(zip(pos.tolist(), w.tolist())) == \
        [(173, 0.5), (173, 1.0), (174, 0.5), (225, 0.5), (225, 1.0), (226, 0.5)]


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


def test_negative_bin(capsys):
    # right read ending at 40, tlen -100: fragment [-60, 40), centre -11, bin -1
    f = F([[100, 150, 147], [20, 40, -100]])
    rows = counts.bincount({'chr1': f}, {'chr1': 2500}, ['chr1'], 1000, tlen=False)
    assert [r[3] for r in rows] == [1, 0, 0] and 'dropped 1 fragment' in capsys.readouterr().err
    rows = counts.bincount({'chr1': f}, {'chr1': 2500}, ['chr1'], 1000, tlen=False, park=True)
    assert [r[3] for r in rows] == [1, 0, 1]                # Park: wraps to the last bin


def test_peaks_greedy():
    # dyad 1000 (3 reads) beats 1050 (2 reads, 50 bp away); 1200 (1 read) is 200 bp away, kept
    rows = [[1000 - 73, 1000, 147]] * 3 + [[1050 - 73, 1050, 147]] * 2 + [[1200 - 73, 1200, 147]]
    pk = counts.peaks(F(rows))
    assert pk == {1000: 3.0, 1200: 1.0}
    assert counts.peaks(F(rows), park=True) == pk


def test_peak_ties():
    # dyads 1000 and 1050 tie (2 reads each) and are 50 bp apart, so only one is kept: by default the
    # lower position; Park mode keeps the one Python 2's sort put first (here the higher one)
    rows = [[1000 - 73, 1000, 147]] * 2 + [[1050 - 73, 1050, 147]] * 2
    assert counts.peaks(F(rows)) == {1000: 2.0}
    assert counts.peaks(F(rows), park=True) == {1050: 2.0}


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
    (header, rows), = score.score([(t, c)], survival_fraction=0.5, missing_as_nan=False)
    # Park: rows scored are control rows with control > 0 (chr1:0, chr1:10); the test-only row is dropped
    T, C = (5 + 1) + (0 + 1), (10 + 1) + (4 + 1)
    want0 = -math.log(((5 + 1) / T) / ((10 + 1) / C)) - math.log(0.5)
    assert [r[1] for r in rows] == ['0', '10']
    assert float(rows[0][3]) == round(want0, 5)
    # revised: the control row without reads stays, as NA; the other values do not change
    (header, rows_na), = score.score([(t, c)], survival_fraction=0.5, missing_as_nan=True)
    assert rows_na[:2] == rows and rows_na[2] == ['chr1', '20', '30', 'NA']


def test_na_rows_downstream(tmp_path):
    # default missing_as_nan: NA rows (bins without input reads) pass through concat and the QC summaries
    import json
    from condense_calc import qc
    hdr, n, rng = ['Chromosome', 'Start', 'End'], 20, np.random.default_rng(1)

    def tab(path, cols):                                   # cols: {name: n values}
        gtab.write(path, hdr + list(cols), [['chr1', 10 * i, 10 * i + 10] + [v[i] for v in cols.values()]
                                            for i in range(n)])
        return path
    score_of = lambda na: ['NA' if i in na else round(float(rng.normal()), 5) for i in range(n)]    # as score() writes
    pts = [tab(tmp_path / f'X_{t}_score.gtab.gz', {f'X_{t}': score_of((3, 7))}) for t in (1, 2)]
    header, rows = score.concatenate(pts)
    gtab.write(tmp_path / 'X.gtab.gz', header, rows)
    rows = [l.split('\t') for l in gtab.open_text(tmp_path / 'X.gtab.gz').read().splitlines()[1:]]
    assert len(rows) == n and rows[3][3:] == rows[7][3:] == ['NA', 'NA']
    tab(tmp_path / 'Y.gtab.gz', {'X_1': score_of((5,)), 'X_2': score_of((5,))})   # a replicate
    bins = [tab(tmp_path / f'X_{t}_bin.gtab.gz', {f'X_{t}': [0 if t == 0 and i in (3, 7) else int(rng.integers(1, 50))
                                                             for i in range(n)]}) for t in (0, 1, 2)]
    samples = [tmp_path / f'X_{t}.json' for t in (0, 1, 2)]
    for t, f in enumerate(samples):
        f.write_text(json.dumps({'sample': f'X_{t}'}))
    (tmp_path / 't.tsv').write_text('tnum\tconc\tsurvival_fraction\n0\t0\t1\n1\t0.2\t0.8\n2\t0.5\t0.4\n')
    d = qc.series_stats('X', samples, bins, tmp_path / 'X.gtab.gz', tmp_path / 't.tsv')
    assert d['zero_input_fraction'] == 2 / n and [p['bins_scored'] for p in d['points']] == [n - 2, n - 2]
    assert all(np.isfinite([p['median'], p['q05'], p['rho_input'], p['rho_supernatant']]).all() for p in d['points'])
    r = qc.replicate_stats('X', {'10kb': [tmp_path / 'X.gtab.gz', tmp_path / 'Y.gtab.gz']})
    assert [p['bins'] for p in r['pairs']] == [n - 3, n - 3] and all(np.isfinite(p['spearman']) for p in r['pairs'])


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


def test_mate_of():
    assert sheet.mate_of('a/x_S1_L001_R1_001.fastq.gz') == 'a/x_S1_L001_R2_001.fastq.gz'
    assert sheet.mate_of('SRR1_1.fastq.gz') == 'SRR1_2.fastq.gz'
    with pytest.raises(sheet.SheetError):
        sheet.mate_of('reads.fastq.gz')


def test_simple_sheet(tmp_path):
    for n in ['in', 'a_L1', 'a_L2', 'b']:
        for r in (1, 2):
            (tmp_path / f'{n}_R{r}.fastq.gz').touch()
    (tmp_path / 's.tsv').write_text('series\tconc\tfraction\tfastq\n'
                                    'X\t0\t\tin_R1.fastq.gz\n'
                                    'X\t0.5\t0.4\ta_L*_R1.fastq.gz\n'
                                    'X\t0.2\t0.9\tb_R1.fastq.gz\n')
    S = sheet.read_simple(tmp_path / 's.tsv', 'hg38', workdir=tmp_path / 'run')
    assert list(S.index) == ['X_0', 'X_1', 'X_2']
    assert S.loc['X_2', 'fastq1'].count(',') == 1 and S.loc['X_2', 'fastq2'].endswith('a_L2_R2.fastq.gz')
    assert (tmp_path / 'run/titration/X.tsv').read_text() == \
        'tnum\tconc\tsurvival_fraction\n0\t0\t1\n1\t0.2\t0.9\n2\t0.5\t0.4\n'
    (tmp_path / 's.tsv').write_text('series\tconc\tfraction\tfastq\nX\t0.5\t1.4\tb_R1.fastq.gz\n')
    with pytest.raises(sheet.SheetError, match='fraction'):
        sheet.read_simple(tmp_path / 's.tsv', 'hg38', workdir=tmp_path / 'run')


def _slurm_cmds(tmp_path, monkeypatch, extra, in_job):
    """The commands `condense-calc run --slurm` would start (nothing is run)."""
    from condense_calc import cli
    cfg = tmp_path / 'config.yaml'
    cfg.write_text('bowtie2: bowtie2\n')                        # this env's: nothing to install
    calls = []
    monkeypatch.setattr(cli.subprocess, 'call', lambda cmd, **k: calls.append(cmd) or 0)
    monkeypatch.setattr(cli.subprocess, 'run',
                        lambda cmd, **k: calls.append(cmd) or type('R', (), {'stdout': '1'})())
    if in_job:                                   # the controller job: runs the workflow itself
        monkeypatch.setenv('SLURM_JOB_ID', '1')
        monkeypatch.setenv('CONDENSE_CONTROLLER', '1')
    else:
        monkeypatch.delenv('SLURM_JOB_ID', raising=False)
        monkeypatch.delenv('CONDENSE_CONTROLLER', raising=False)
    try:
        cli.main(['run', str(cfg), '-o', str(tmp_path / 'out'), '--slurm', *extra])
    except SystemExit:
        pass
    return calls


def test_slurm_run_keeps_rule_threads(tmp_path, monkeypatch):
    # inside the 1-core controller job, Snakemake must get a thread budget (--cores); without it,
    # every rule's threads are capped at the controller's 1 CPU (seen 2026-10-18: align on 1 CPU)
    cmd = _slurm_cmds(tmp_path, monkeypatch, [], in_job=True)[0]
    assert '--executor' in cmd and '--cores' in cmd
    assert int(cmd[cmd.index('--cores') + 1]) >= 64


def test_run_finds_tools_of_its_environment(tmp_path, monkeypatch):
    # started as /path/env/bin/python -m condense_calc without activating the env: snakemake must still be found
    import sys
    from pathlib import Path
    from condense_calc import cli
    cfg = tmp_path / 'config.yaml'
    cfg.write_text('bowtie2: bowtie2\n')                        # this env's: nothing to install
    monkeypatch.setenv('SLURM_JOB_ID', '1')
    monkeypatch.setenv('CONDENSE_CONTROLLER', '1')
    monkeypatch.setenv('PATH', '/usr/bin')
    envs = []
    monkeypatch.setattr(cli.subprocess, 'call', lambda cmd, env=None, **k: envs.append(env) or 0)
    with pytest.raises(SystemExit):
        cli.main(['run', str(cfg), '-o', str(tmp_path / 'out'), '--slurm'])
    assert envs[0]['PATH'].split(':')[0] == str(Path(sys.executable).parent)


def test_slurm_controller_mail(tmp_path, monkeypatch):
    cmd = _slurm_cmds(tmp_path, monkeypatch, ['--mail', 'me@example.org'], in_job=False)[0]
    assert cmd[0] == 'sbatch' and '--mail-user=me@example.org' in cmd
    assert '--mail-type=FAIL,TIME_LIMIT_80' in cmd
    assert not any(c.startswith('--mail') for c in _slurm_cmds(tmp_path, monkeypatch, [], in_job=False)[0])


def _fai(tmp_path, seqs):
    p = tmp_path / 'g.fa.fai'
    p.write_text(''.join(f'{n}\t{L}\t0\t60\t61\n' for n, L in seqs))
    return p


def test_default_chroms(tmp_path):
    from condense_calc.reference import default_chroms
    hg = [(f'chr{c}', 50_000_000) for c in list(range(1, 23)) + ['X', 'Y']]
    hg += [('chrM', 16_569), ('chrUn_KI270302v1', 2_274), ('chrEBV', 171_823)]
    assert default_chroms(_fai(tmp_path, hg)) == [n for n, _ in hg[:24]]          # as the GEO run
    ucsc = sorted(hg[:24]) + hg[24:]                                              # UCSC hg38.fa: chr1, chr10, ...
    assert default_chroms(_fai(tmp_path, ucsc)) == [n for n, _ in hg[:24]]
    dm6 = [('chr2L', 23_513_712), ('chr2R', 25_286_936), ('chr3L', 28_110_227), ('chr3R', 32_079_331),
           ('chr4', 1_348_131), ('chrX', 23_542_271), ('chrY', 3_667_352), ('chrM', 19_524),
           ('chrUn_DS483562v1', 1_001)]
    assert default_chroms(_fai(tmp_path, dm6)) == [n for n, _ in dm6[:7]]
    yeast = [(f'chr{r}', 800_000) for r in 'I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI'.split()]
    assert default_chroms(_fai(tmp_path, yeast + [('chrM', 85_779)])) == [n for n, _ in yeast]
    ncbi = [('NC_000001.11', 248_956_422), ('NC_000002.12', 242_193_529), ('NC_012920.1', 16_569)]
    assert default_chroms(_fai(tmp_path, ncbi)) == ['NC_000001.11', 'NC_000002.12']
    with pytest.raises(ValueError):
        default_chroms(_fai(tmp_path, [('chrM', 16_569)]))


def test_slurm_from_interactive_job_submits_controller(tmp_path, monkeypatch):
    # inside salloc / Open OnDemand (SLURM_JOB_ID set, not the controller): submit, do not run in the session
    monkeypatch.setenv('SLURM_JOB_ID', '7')
    from condense_calc import cli
    cfg = tmp_path / 'config.yaml'
    cfg.write_text('bowtie2: bowtie2\n')                        # this env's: nothing to install
    seen = []
    monkeypatch.setattr(cli.subprocess, 'call', lambda cmd, **k: seen.append(cmd) or 0)
    monkeypatch.setattr(cli.subprocess, 'run', lambda cmd, **k: seen.append(cmd) or type('R', (), {'stdout': '1'})())
    monkeypatch.delenv('CONDENSE_CONTROLLER', raising=False)
    cli.main(['run', str(cfg), '-o', str(tmp_path / 'out'), '--slurm'])
    assert seen[0][0] == 'sbatch' and '--export=ALL,CONDENSE_CONTROLLER=1' in seen[0]


def _local_cmd(tmp_path, monkeypatch, extra, env):
    from condense_calc import cli
    cfg = tmp_path / 'config.yaml'
    cfg.write_text('bowtie2: bowtie2\n')                        # this env's: nothing to install
    for k in ('SLURM_JOB_ID', 'SLURM_MEM_PER_NODE', 'SLURM_MEM_PER_CPU', 'CONDENSE_CONTROLLER'):
        monkeypatch.delenv(k, raising=False)
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    seen = []
    monkeypatch.setattr(cli.subprocess, 'call', lambda cmd, **k: seen.append(cmd) or 0)
    try:
        cli.main(['run', str(cfg), '-o', str(tmp_path / 'out'), *extra])
    except SystemExit as e:
        if not seen:
            return str(e.code)
    return seen[0]


def test_local_run_resources(tmp_path, monkeypatch):
    from condense_calc import cli
    monkeypatch.setattr(cli.os, 'sched_getaffinity', lambda pid: {0, 1, 2})
    cmd = _local_cmd(tmp_path, monkeypatch, [], {'SLURM_JOB_ID': '5', 'SLURM_MEM_PER_NODE': '4000'})
    assert cmd[cmd.index('--cores') + 1] == '3'                   # the job's CPUs, not the node's
    assert cmd[cmd.index('--resources') + 1] == 'mem_mb=4000'
    cmd = _local_cmd(tmp_path, monkeypatch, ['--cores', '2'],
                     {'SLURM_JOB_ID': '5', 'SLURM_MEM_PER_CPU': '1500'})
    assert cmd[cmd.index('--cores') + 1] == '2' and cmd[cmd.index('--resources') + 1] == 'mem_mb=3000'


def test_local_run_on_login_node_refused(tmp_path, monkeypatch):
    from condense_calc import cli
    monkeypatch.setattr(cli.shutil, 'which', lambda name: '/usr/bin/sbatch')
    msg = _local_cmd(tmp_path, monkeypatch, [], {})
    assert isinstance(msg, str) and 'login node' in msg
    cmd = _local_cmd(tmp_path, monkeypatch, ['--cores', '2'], {})          # explicit: allowed
    assert cmd[cmd.index('--cores') + 1] == '2'


def test_dry_run_does_not_fetch_genome(tmp_path, monkeypatch, capsys):
    from condense_calc import cli, reference
    sh = tmp_path / 'sheet.tsv'
    sh.write_text('series\tconc\tfraction\tfastq\n')
    monkeypatch.setattr(cli, '_check', lambda *a, **k: None)
    monkeypatch.setattr(reference, 'missing', lambda g, r=None: 'download hg38 into X')
    monkeypatch.setattr(reference, 'ensure', lambda *a, **k: (_ for _ in ()).throw(AssertionError('fetched')))
    with pytest.raises(SystemExit) as e:
        cli.main(['run', str(sh), '-o', str(tmp_path / 'out'), '--slurm', '-n'])
    assert e.value.code == 0 and 'would first download hg38' in capsys.readouterr().out


def test_genome_prepared_in_own_job(tmp_path, monkeypatch):
    from condense_calc import cli, reference
    sh = tmp_path / 'sheet.tsv'
    sh.write_text('series\tconc\tfraction\tfastq\n')
    for k in ('SLURM_JOB_ID', 'CONDENSE_CONTROLLER'):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setattr(cli, '_check', lambda *a, **k: None)
    monkeypatch.setattr(reference, 'missing', lambda g, r=None: 'download hg38 into X')
    ids = iter(['111', '222'])
    seen = []
    monkeypatch.setattr(cli.subprocess, 'run',
                        lambda cmd, **k: seen.append(cmd) or type('R', (), {'stdout': next(ids)})())
    cli.main(['run', str(sh), '-o', str(tmp_path / 'out'), '--slurm', '--revised'])
    ref, ctl = seen
    opt = lambda cmd, o: cmd[cmd.index(o) + 1]
    assert 'reference' in ref[-1] and opt(ref, '-c') == '1' and opt(ref, '-t') == '0-12:00'    # hg38: a download
    assert '--dependency=afterok:111' in ctl and opt(ctl, '-c') == '1'
    # an index to build (the default hg38_ucsc, or a FASTA): 8 cores
    seen.clear()
    ids = iter(['111', '222'])
    cli.main(['run', str(sh), '-o', str(tmp_path / 'out'), '--slurm'])
    assert opt(seen[0], '-c') == '8' and opt(seen[0], '--mem') == '32G' and '--threads 8' in seen[0][-1]
    # and not on a login node by hand
    monkeypatch.setattr(cli, '_on_login_node', lambda: True)
    with pytest.raises(SystemExit, match='not light work'):
        cli.main(['reference', 'hg38_ucsc'])


def test_settings_resolve_and_validate():
    from condense_calc import settings
    S = settings.resolve({})
    assert S == settings.DEFAULTS                                            # Park et al.'s processing
    S = settings.resolve({'reads': {'dedup': 'true'}, 'align_threads': 24, 'samples': 'x.tsv', 'qc': 'OFF'})
    assert S['reads'] == {**settings.DEFAULTS['reads'], 'dedup': True}      # merged one level deep
    assert S['threads'] == 24 and S['qc'] == 'off'                          # earlier name of threads
    assert settings.resolve({'park': 'false'})['park'] is False              # as Snakemake's --config park=false
    with pytest.raises(settings.SettingsError, match="did you mean 'reads'"):
        settings.resolve({'read': {'min_mapq': 10}})
    with pytest.raises(settings.SettingsError, match="reads.min_mapqq"):
        settings.resolve({'reads': {'min_mapqq': 10}})
    with pytest.raises(settings.SettingsError, match='qc_thresholds.min_replicate_rh'):
        settings.resolve({'qc_thresholds': {'min_replicate_rh': 0.7}})
    with pytest.raises(settings.SettingsError, match='qc must be one of'):
        settings.resolve({'qc': 'full'})
    import yaml
    S = settings.resolve(yaml.safe_load('fastq_check: off\nqc: off\n'))     # YAML reads a bare off as false
    assert S['fastq_check'] == 'off' and S['qc'] == 'off'
    bad = [({'reads': {'min_mapq': 'ten'}}, 'reads.min_mapq must be an integer'),
           ({'threads': 'eight'}, 'threads must be an integer'), ({'max_runtime': '1d'}, 'max_runtime must be an'),
           ({'fit': {'min_rsq': 'high'}}, 'fit.min_rsq must be a number'), ({'mscale': {'Ca': 'x'}}, 'mscale.Ca'),
           ({'bin_sizes': 1000}, 'bin_sizes must be a list'), ({'bin_sizes': [0]}, 'bin_sizes must be positive'),
           ({'fit': {'model': 'logistic'}}, 'fit.model must be one of'),
           ({'cells': {'GM': {'genome': 'hg38', 'chrom': ['chr1']}}}, "did you mean 'cells.GM.chroms'")]
    for cfg, msg in bad:
        with pytest.raises(settings.SettingsError, match=msg):
            settings.resolve(cfg)
    assert settings.resolve({'reads': {'min_mapq': '10'}})['reads']['min_mapq'] == 10    # --config gives text


def test_stage_defaults_from_settings():
    import argparse
    from condense_calc import cli, settings
    a = argparse.Namespace(park=True, mm=None, min_mapq=None, dedup=None, missing_as_nan=True)
    cli._stage_defaults(a)
    R = settings.DEFAULTS['reads']
    assert (a.mm, a.min_mapq, a.dedup) == (R['mm'], R['min_mapq'], R['dedup'])
    assert a.missing_as_nan is True                                          # given: kept


def test_run_stops_on_unknown_setting(tmp_path, monkeypatch):
    from condense_calc import cli
    cfg = tmp_path / 'config.yaml'
    cfg.write_text('reads: {min_mapq: 10}\nbin_size: [1000]\n')
    monkeypatch.setattr(cli.subprocess, 'call', lambda *a, **k: pytest.fail('started'))
    with pytest.raises(SystemExit) as e:
        cli.main(['run', str(cfg), '-o', str(tmp_path / 'out'), '--cores', '1'])
    assert "unknown setting 'bin_size'; did you mean 'bin_sizes'?" in str(e.value.code)
    assert not (tmp_path / 'out').exists()


def test_mode_reaches_the_workflow(tmp_path, monkeypatch):
    # a config run gets --config park=false (--revised) or park=true (--park); a sheet run has it in condense.yaml
    # (test_park_maps_hg38)
    from condense_calc import cli, reference
    monkeypatch.setattr(reference, 'missing_bowtie2', lambda *a, **k: None)   # bowtie2 2.4.1 is there
    monkeypatch.setattr(reference, 'ensure_bowtie2', lambda *a, **k: None)
    for flag, conf in (('--revised', 'park=false'), ('--park', 'park=true')):
        cmd = _local_cmd(tmp_path, monkeypatch, ['--cores', '1', flag], {})
        assert cmd[cmd.index('--config') + 1:cmd.index('--config') + 2] == [conf]
    assert not any(c.startswith('park=') for c in _local_cmd(tmp_path, monkeypatch, ['--cores', '1'], {}))
    with pytest.raises(SystemExit):
        cli.main(['run', 'x.tsv', '--park', '--revised'])                    # one or the other


def test_config_refs_used_for_bowtie2(tmp_path, monkeypatch):
    # a config's own `refs:` (relative to the config file) is where the CLI looks for bowtie2, as the workflow does
    from condense_calc import cli, reference
    cfg = tmp_path / 'config.yaml'
    cfg.write_text('park: true\nrefs: myrefs\n')
    seen = []
    monkeypatch.setattr(reference, 'missing_bowtie2', lambda s, r=None: seen.append(r))
    monkeypatch.setattr(reference, 'ensure_bowtie2', lambda s, r=None: seen.append(r))
    monkeypatch.setattr(cli.subprocess, 'call', lambda cmd, **k: seen.append(cmd) or 0)
    with pytest.raises(SystemExit):
        cli.main(['run', str(cfg), '-o', str(tmp_path / 'out'), '--cores', '1'])
    assert seen[:-1] == [str(tmp_path / 'myrefs')] * (len(seen) - 1) and f'refs={tmp_path}/myrefs' in seen[-1]


def _sheet_config(tmp_path, monkeypatch, extra, settings_text=None):
    """condense.yaml as `condense SHEET` writes it (nothing is fetched or run)."""
    import yaml
    from condense_calc import cli, reference
    sh = tmp_path / 'sheet.tsv'
    sh.write_text('series\tconc\tfraction\tfastq\n')
    if settings_text is not None:
        (tmp_path / 'settings.yaml').write_text(settings_text)
        extra = extra + ['--settings', str(tmp_path / 'settings.yaml')]
    monkeypatch.setattr(cli, '_check', lambda *a, **k: None)
    monkeypatch.setattr(reference, 'ensure', lambda g, *a, **k: (g, f'/refs/{g}/{g}.fa', f'/refs/{g}/{g}'))
    monkeypatch.setattr(reference, 'missing_bowtie2', lambda *a, **k: None)   # bowtie2 2.4.1 is there
    monkeypatch.setattr(reference, 'ensure_bowtie2', lambda *a, **k: None)
    monkeypatch.setattr(cli.subprocess, 'call', lambda cmd, **k: 0)
    with pytest.raises(SystemExit):
        cli.main(['run', str(sh), '-o', str(tmp_path / 'out'), '--cores', '1', *extra])
    return yaml.safe_load((tmp_path / 'out' / 'condense.yaml').read_text())


def test_settings_file_nucleosomes_kept(tmp_path, monkeypatch):
    # nucleosomes: true in --settings was overwritten by the flag's default (False)
    assert _sheet_config(tmp_path, monkeypatch, [], 'nucleosomes: true\n')['nucleosomes'] is True
    assert 'nucleosomes' not in _sheet_config(tmp_path, monkeypatch, [], 'threads: 2\n')
    assert _sheet_config(tmp_path, monkeypatch, ['--nucleosomes'])['nucleosomes'] is True


def test_sheet_titration_verbatim(tmp_path):
    # the sheet's own strings: {f:g} wrote 6 digits (0.097744361 -> 0.0977444), which cost GEO byte identity
    for n in ['in', 'a']:
        for r in (1, 2):
            (tmp_path / f'{n}_R{r}.fastq.gz').touch()
    (tmp_path / 's.tsv').write_text('series\tconc\tfraction\tfastq\n'
                                    'X\t0\t\tin_R1.fastq.gz\n'
                                    'X\t0.79432823\t0.097744361\ta_R1.fastq.gz\n')
    notes = []
    sheet.read_simple(tmp_path / 's.tsv', 'hg38', workdir=tmp_path / 'run', note=notes.append)
    assert (tmp_path / 'run/titration/X.tsv').read_text() == \
        'tnum\tconc\tsurvival_fraction\n0\t0\t1\n1\t0.79432823\t0.097744361\n'
    assert not notes
    (tmp_path / 's.tsv').write_text('series\tconc\tfraction\tfastq\nX\t0\t1\tin_R1.fastq.gz\n'
                                    'X\t0.794\t0.0977444\ta_R1.fastq.gz\n')
    sheet.read_simple(tmp_path / 's.tsv', 'hg38', workdir=tmp_path / 'run', note=notes.append)
    assert len(notes) == 1 and '0.0977444' in notes[0]
    assert sheet._sig_digits('0.400') == 3 and sheet._sig_digits('9.7744361e-2') == 8


def test_revised_preset_and_explicit_keys():
    import argparse
    from condense_calc import cli, settings
    assert settings.resolve({})['reads']['min_mapq'] == 0
    rev = settings.resolve({'park': False})
    assert rev['reads'] == {'mm': 10, 'min_mapq': 10, 'dedup': False}
    assert (rev['score'], rev['bowtie2'], rev['fastq_check']) == ({'missing_as_nan': True}, 'bowtie2', 'lite')
    assert settings.resolve({'park': False, 'reads': {'min_mapq': 5}})['reads']['min_mapq'] == 5    # explicit wins
    a = argparse.Namespace(park=False, min_mapq=None, mm=None)
    cli._stage_defaults(a)
    assert (a.min_mapq, a.mm) == (10, 10)
    assert settings.describe('reads', 'min_mapq') == 'default 0; --revised: 10'


def test_bdata_opt_in():
    from condense_calc import settings
    assert settings.resolve({})['bdata'] is False and settings.resolve({'park': False})['bdata'] is False
    assert settings.resolve({'bdata': 'true'})['bdata'] is True


def test_reference_presets(tmp_path):
    import hashlib
    from condense_calc import reference
    fai = _fai(tmp_path, [('chr1', 1000), ('chr2', 500)])
    assert reference.fingerprint(fai) == hashlib.md5(b'chr1\t1000\nchr2\t500\n').hexdigest()   # cut -f1,2 | md5sum
    assert reference.locate('hg38_ucsc', tmp_path)[1:] == (tmp_path / 'hg38_ucsc/hg38_ucsc.fa', tmp_path / 'hg38_ucsc/hg38_ucsc')
    assert reference.missing('mm10', tmp_path).startswith('download mm10')
    # an older install under a preset's name (here 2 sequences) is refused, not used
    d = tmp_path / 'mm10'
    d.mkdir()
    (d / 'mm10.fa.fai').write_text(fai.read_text())
    for f in (reference.missing, reference.ensure):
        with pytest.raises(SystemExit, match='(?s)holds another reference .*mv '):
            f('mm10', tmp_path)


def test_park_maps_hg38(tmp_path, monkeypatch):
    cfg = _sheet_config(tmp_path, monkeypatch, [])
    assert 'park' not in cfg and cfg['genome'] == 'hg38_ucsc'      # the default: Park et al.'s hg38
    assert _sheet_config(tmp_path, monkeypatch, ['--park'])['genome'] == 'hg38_ucsc'
    cfg = _sheet_config(tmp_path, monkeypatch, ['-g', 'mm10'])
    assert cfg['genome'] == 'mm10'                                  # Park's mouse reference is the default one
    cfg = _sheet_config(tmp_path, monkeypatch, ['--revised'])
    assert cfg['park'] is False and cfg['genome'] == 'hg38'
    assert _sheet_config(tmp_path, monkeypatch, [], 'park: false\n')['genome'] == 'hg38'


def test_bowtie2_spec(tmp_path, monkeypatch):
    import yaml
    from condense_calc import reference, settings
    assert reference.bowtie2_dir('bowtie2') is None and reference.missing_bowtie2('bowtie2', tmp_path) is None
    env = tmp_path / 'envs/bowtie2-2.4.1'
    assert reference.bowtie2_dir('2.4.1', tmp_path) == str(env / 'bin')
    assert reference.missing_bowtie2('2.4.1', tmp_path).startswith('install bowtie2 2.4.1')
    # install (subprocess mocked: nothing is installed); bioconda's 2.4.1 prints no version, conda-meta has it
    calls = []

    def run(cmd, **k):
        calls.append(cmd)
        if cmd[1] == 'create':
            (env / 'bin').mkdir(parents=True)
            (env / 'bin/bowtie2').touch()
            (env / 'conda-meta').mkdir()
            (env / 'conda-meta/bowtie2-2.4.1-py38h1c8e9b9_3.json').write_text('{"version": "2.4.1"}')
        return type('R', (), {'stdout': f'{env}/bin/bowtie2-align-s version \n64-bit\n'})()
    monkeypatch.setattr(reference.subprocess, 'run', run)
    monkeypatch.setenv('MAMBA_EXE', '/x/mamba')
    assert reference.ensure_bowtie2('2.4.1', tmp_path) == str(env / 'bin')
    assert calls[0] == ['/x/mamba', 'create', '-y', '-q', '-p', str(env), '-c', 'conda-forge', '-c', 'bioconda',
                        'bowtie2=2.4.1']
    assert calls[1] == ['bowtie2', '--version'] and len(calls) == 2
    assert reference.missing_bowtie2('2.4.1', tmp_path) is None
    assert reference.bowtie2_dir(str(env)) == reference.bowtie2_dir(str(env / 'bin/bowtie2')) == str(env / 'bin')
    with pytest.raises(SystemExit, match='no bowtie2 in'):
        reference.missing_bowtie2(str(tmp_path / 'nothing/bowtie2'))
    # a partial version would install 2.5.5 and then never match '2.5'
    with pytest.raises(SystemExit, match='full version'):
        reference.missing_bowtie2('2.5', tmp_path)
    with pytest.raises(settings.SettingsError, match='full version'):
        settings.resolve(yaml.safe_load('bowtie2: 2.10'))
    assert settings.resolve(yaml.safe_load('bowtie2: 2.4.1'))['bowtie2'] == '2.4.1'


def test_park_bowtie2_in_reference_job(tmp_path, monkeypatch):
    # Park mode (the default): bowtie2 2.4.1 is installed by the genome-preparation job, not on the login node or the controller
    from condense_calc import cli, reference
    sh = tmp_path / 'sheet.tsv'
    sh.write_text('series\tconc\tfraction\tfastq\n')
    for k in ('SLURM_JOB_ID', 'CONDENSE_CONTROLLER'):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setattr(cli, '_check', lambda *a, **k: None)
    monkeypatch.setattr(reference, 'missing', lambda g, r=None: None)
    monkeypatch.setattr(reference, 'missing_bowtie2', lambda s, r=None: None if s == 'bowtie2' else f'install {s}')
    ids = iter(['111', '222'])
    seen = []
    monkeypatch.setattr(cli.subprocess, 'run',
                        lambda cmd, **k: seen.append(cmd) or type('R', (), {'stdout': next(ids)})())
    cli.main(['run', str(sh), '-o', str(tmp_path / 'out'), '--slurm'])
    ref, ctl = seen
    assert '--bowtie2 2.4.1' in ref[-1] and ref[-1].endswith('hg38_ucsc') and '--dependency=afterok:111' in ctl
    assert ref[ref.index('-c') + 1] == '1' and ref[ref.index('-t') + 1] == '0-02:00'    # only a conda install
    # the controller itself never installs it
    monkeypatch.setenv('SLURM_JOB_ID', '222')
    monkeypatch.setenv('CONDENSE_CONTROLLER', '1')
    with pytest.raises(SystemExit, match='First install 2.4.1'):
        cli.main(['run', str(sh), '-o', str(tmp_path / 'out'), '--slurm', '--park'])


def _fq(path, names, qual='I' * 8, gz=False):
    import gzip
    text = ''.join(f'{n}\nACGTACGT\n+\n{qual}\n' for n in names)
    (gzip.open(path, 'wt') if gz else open(path, 'w')).write(text)
    return path


def test_fastq_checks(tmp_path):
    from condense_calc import fastq
    inst = _fq(tmp_path / 'i_1.fastq.gz', ['@VL1:75:AAG:1:1101:1:1 1:N:0:ATTACTCG+AGGTTAGG'] * 3, gz=True)
    bare = _fq(tmp_path / 'b_2.fastq', ['@SRR1.1 2', '@SRR1.2 2'])                      # '@<sn> <mate>': fine
    ena = _fq(tmp_path / 'e_1.fastq.gz', ['@SRR31627121.1 VL1:75:AAG:1:1101:1:1/1'], gz=True)
    two = _fq(tmp_path / 't_1.fastq', ['@VL1:75 1:N:0:A', '@VL1:76  1:N:0:A'])           # two spaces, read 2
    lite = _fq(tmp_path / 'l_1.fastq.gz', ['@x 1:N:0:A'] * 4, qual='?' * 8, gz=True)
    ok = lambda f, mate, mode: not fastq.check_file(f, mate, mode)['problems']
    assert ok(inst, 1, 'park') and ok(bare, 2, 'park') and not ok(inst, 2, 'park')         # mate number checked
    assert ok(ena, 1, 'lite') and not ok(ena, 1, 'park') and not ok(two, 1, 'park')
    assert not ok(lite, 1, 'lite') and ok(lite, 1, 'off')
    r = fastq.check_file(lite, 1, 'lite')
    assert r['sra_lite'] and r['reads_checked'] == 4 and 'condense fetch' in r['problems'][0]
    assert fastq.check([inst], [bare], 'off') == []
    assert len(fastq.problems(fastq.check([ena], [bare], 'park'))) == 1
    # SRA Lite and ENA names at once: both messages give the setting that aligns the file anyway, off
    both = _fq(tmp_path / 'x_1.fastq.gz', ['@SRR1.1 VL1:75:AAG:1:1101:1:1/1'], qual='?' * 8, gz=True)
    msgs = fastq.check_file(both, 1, 'park')['problems']
    assert len(msgs) == 2 and all('set fastq_check: off (in the --settings file' in m for m in msgs)
    assert 'set fastq_check: lite' in fastq.check_file(ena, 1, 'park')['problems'][0]
    (tmp_path / 'n.fq').write_text('@a 1\nAC\n+\n??\n' * 2 + '@a 1\nAC\n+\nII\n')     # only n reads are read
    assert fastq.check_file(tmp_path / 'n.fq', 1, 'lite', n=2)['sra_lite']


def test_sheet_run_stops_on_sra_lite(tmp_path, monkeypatch):
    from condense_calc import cli
    for n in ['in', 'a']:
        for m in (1, 2):
            _fq(tmp_path / f'{n}_R{m}.fastq.gz', [f'@x {m}:N:0:A'] * 3, qual='?' * 8 if n == 'a' else 'I' * 8, gz=True)
    (tmp_path / 's.tsv').write_text('series\tconc\tfraction\tfastq\nX\t0\t1\tin_R1.fastq.gz\nX\t0.5\t0.4\ta_R1.fastq.gz\n')
    monkeypatch.setattr(cli.subprocess, 'run', lambda *a, **k: pytest.fail('submitted'))
    with pytest.raises(SystemExit) as e:
        cli.main(['run', str(tmp_path / 's.tsv'), '-o', str(tmp_path / 'out'), '--slurm'])
    assert 'a_R1.fastq.gz' in str(e.value.code) and 'SRA Lite' in str(e.value.code)
    (tmp_path / 'out/bam').mkdir()
    (tmp_path / 'out/bam/X_1.bam').touch()               # already aligned: not checked again
    cli._check(tmp_path / 's.tsv', 'lite', tmp_path / 'out')


def test_config_samples_errors(tmp_path):
    # a wrong samples path or a missing column: one clear line, not a pandas traceback
    from condense_calc import cli
    cfg = tmp_path / 'c.yaml'
    cfg.write_text('samples: nosuch.tsv\n')
    with pytest.raises(sheet.SheetError, match='samples nosuch.tsv: .*No such file'):
        cli._check(cfg, 'lite')
    (tmp_path / 's.tsv').write_text('cell\tsample\tagent\ttnum\trep\tdepth\tfastq1\nH1\tNCP\tsp\t0\t1\t\ta.fq\n')
    cfg.write_text('samples: s.tsv\n')
    with pytest.raises(sheet.SheetError, match="samples s.tsv: no column 'fastq2'"):
        cli._check(cfg, 'lite')


def test_config_run_rejects_sheet_options(tmp_path):
    # --settings, -g and --nucleosomes would be ignored by a config run: stop instead
    from condense_calc import cli
    cfg = tmp_path / 'c.yaml'
    cfg.write_text('samples: s.tsv\n')
    for extra in (['--settings', 'x.yaml'], ['-g', 'mm10'], ['--nucleosomes']):
        with pytest.raises(SystemExit, match='not used in a config run'):
            cli.main(['run', str(cfg), '-o', str(tmp_path / 'out'), *extra])


def test_park_chry_rule():
    from condense_calc import settings
    hs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    for label in ('GM', 'GM_NCP_sp', 'mCD8T', 'mCD8T:WT', 'mCD8T_DFMO_NCP_sp', 'mCD8T-ODCKO'):
        assert settings.park_chroms(label, hs) == hs[:-1], label
    for label in ('H1', 'E14_NCP_sp', 'GMX', 'H1_GM'):
        assert settings.park_chroms(label, hs) == hs, label
    # a cell's chromosomes: explicit chroms always win; Park mode drops chrY for GM by default only
    default = lambda: hs
    assert settings.cell_chroms('GM', None, default) == (hs, None)
    assert settings.cell_chroms('GM', None, default, park=True) == (hs[:-1], None)
    assert settings.cell_chroms('GM', hs, default, park=True) == (hs, None)
    assert settings.cell_chroms('H1', None, default, park=True) == (hs, None)
    chroms, warning = settings.cell_chroms('GM12878_NCP_sp', None, default, park=True)     # not Park's cell name
    assert chroms == hs and 'GM12878_NCP_sp keeps chrY' in warning
    assert settings.cell_chroms('GM12878_NCP_sp', None, default) == (hs, None)


def test_run_record(tmp_path, monkeypatch):
    import json
    from condense_calc import record, reference, settings
    _fai(tmp_path, [('chr1', 1000)])
    genomes = {'g': {'fasta': str(tmp_path / 'g.fa'), 'bowtie2_index': str(tmp_path / 'g')}}
    logs = tmp_path / 'logs'
    logs.mkdir()
    (logs / 'A_0.bowtie2_version').write_text('/x/bowtie2-align-s version 2.5.5\n64-bit\n')
    (logs / 'A_0.fastq_check.json').write_text('{"mode": "park", "files": []}')
    env = tmp_path / 'envs/bowtie2-2.4.1'                    # bioconda 2.4.1 prints no version; conda-meta has it
    (env / 'conda-meta').mkdir(parents=True)
    (env / 'conda-meta/bowtie2-2.4.1-py38h1c8e9b9_3.json').write_text('{"version": "2.4.1"}')
    (logs / 'A_1.bowtie2_version').write_text(f'{env}/bin/bowtie2-align-s version \n64-bit\n')
    monkeypatch.setattr(reference, 'bowtie2_version', lambda d=None: '2.5.5')
    S = settings.resolve({'park': True, 'bowtie2': 'bowtie2'})
    record.write(tmp_path / 'run_info.json', S, genomes, None, {'A': ['chr1']}, ['A_0', 'A_1'], logs=logs)
    r = json.loads((tmp_path / 'run_info.json').read_text())
    assert r['mode'] == 'park' and r['settings']['reads']['min_mapq'] == 0 and r['bowtie2']['version'] == '2.5.5'
    assert r['references']['g'] == {**genomes['g'], 'sequences': 1, 'preset': None,
                                    'fingerprint': reference.fingerprint(tmp_path / 'g.fa.fai')}
    assert r['libraries']['A_0']['bowtie2'].endswith('2.5.5') and r['libraries']['A_1']['fastq_check'] is None
    assert [r['libraries'][n]['bowtie2_version'] for n in ('A_0', 'A_1')] == ['2.5.5', '2.4.1']
    assert len(r['warnings']) == 1 and 'not a reference Park et al. used' in r['warnings'][0]
    assert record.park_warnings(settings.resolve({'park': False}), genomes) == []
    assert record.code_version().startswith('0.')


def test_fetch_commands(tmp_path, monkeypatch):
    # condense fetch, with prefetch / fasterq-dump / pigz mocked (nothing is downloaded)
    import gzip
    from pathlib import Path
    from condense_calc import cli, fetch
    calls, prefetch_rc, spots = [], iter([1, 0]), {'n': 3}
    R = lambda rc=0, out='': type('R', (), {'returncode': rc, 'stdout': '', 'stderr': out})()

    def run(cmd, **k):
        calls.append(cmd)
        tool, tmp = Path(cmd[0]).name, tmp_path / f'tmp_{cmd[-1].split("/")[-1].split(".")[0]}'
        if tool == 'prefetch':
            rc = next(prefetch_rc, 0)
            if rc == 0:
                (tmp / cmd[-1]).mkdir(parents=True, exist_ok=True)
                (tmp / cmd[-1] / f'{cmd[-1]}.sra').touch()
            return R(rc)
        if tool == 'fasterq-dump':
            run_id = Path(cmd[-1]).stem
            for m in (1, 2):
                _fq(Path(cmd[cmd.index('-O') + 1]) / f'{run_id}_{m}.fastq',
                    [f'@VL1:1:F:1:1:1:{i} {m}:N:0:ACGT' for i in range(3)])
            return R(0, f'spots read      : {spots["n"]}\nreads read      : 6\n')
        f = Path(cmd[-1])                                               # pigz
        if Path(f'{f}.gz').exists():                                    # as pigz and gzip: no overwriting
            raise fetch.subprocess.CalledProcessError(25, cmd)
        gzip.open(f'{f}.gz', 'wb').write(f.read_bytes())
        f.unlink()
        return R()
    monkeypatch.setattr(fetch.subprocess, 'run', run)
    monkeypatch.setattr(fetch.shutil, 'which', lambda n: f'/sw/bin/{n}')
    waits = []
    out = fetch.fetch('SRR1', tmp_path, threads=4, sleep=waits.append)
    assert [p.name for p in out] == ['SRR1_1.fastq.gz', 'SRR1_2.fastq.gz'] and all(p.exists() for p in out)
    assert calls[0] == ['/sw/bin/prefetch', '--verify', 'yes', '--max-size', '200G', '-O', str(tmp_path / 'tmp_SRR1'),
                        'SRR1']
    assert calls[1] == calls[0] and waits == [120]                      # failed once: retried after 2 min
    dump = calls[2]
    assert dump[0] == '/sw/bin/fasterq-dump' and '--split-files' in dump and dump[-1].endswith('SRR1/SRR1.sra')
    assert dump[dump.index('--seq-defline') + 1] == '@$sn $ri:N:0:$sg' and dump[dump.index('--qual-defline') + 1] == '+'
    assert calls[3] == ['/sw/bin/pigz', '-p', '4', str(tmp_path / 'tmp_SRR1/SRR1_1.fastq')]
    assert not (tmp_path / 'tmp_SRR1').exists()                          # .sra and temporary files deleted
    assert fetch.fetch('SRR1', tmp_path) == out and len(calls) == 5      # already there: nothing run
    # a job stopped while compressing left tmp_SRR3/SRR3_1.fastq.gz: the next fetch replaces it
    (tmp_path / 'tmp_SRR3').mkdir()
    (tmp_path / 'tmp_SRR3/SRR3_1.fastq.gz').write_bytes(b'partial')
    assert all(p.exists() for p in fetch.fetch('SRR3', tmp_path))
    # files of the same names from ENA are not taken for ours
    for m in (1, 2):
        _fq(tmp_path / f'SRR4_{m}.fastq.gz', [f'@SRR4.1 VL1:1:F:1:1:1:1/{m}'], gz=True)
    with pytest.raises(fetch.FetchError, match='Move them away'):
        fetch.fetch('SRR4', tmp_path)
    spots['n'] = 4                                                      # fewer reads than the run's spots
    with pytest.raises(fetch.FetchError, match=r'\[3, 3\] reads .* 4 spots'):
        fetch.fetch('SRR2', tmp_path)
    monkeypatch.setattr(cli, '_on_login_node', lambda: False)
    with pytest.raises(SystemExit, match='1 of 2 runs failed: SRR2'):
        cli.main(['fetch', 'SRR1', 'SRR2', '-o', str(tmp_path)])

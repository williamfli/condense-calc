"""Small hand-checkable tests of each stage. Run: pytest tests/test_units.py (in the condense-calc env)."""
import math

import numpy as np
import pytest

from condense_calc import counts, gtab, score, sheet
from condense_calc.fragments import Fragments


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
    (header, rows), = score.score([(t, c)], survival_fraction=0.5, missing_as_nan=False)
    # Park: rows scored are control rows with control > 0 (chr1:0, chr1:10); the test-only row is dropped
    T, C = (5 + 1) + (0 + 1), (10 + 1) + (4 + 1)
    want0 = -math.log(((5 + 1) / T) / ((10 + 1) / C)) - math.log(0.5)
    assert [r[1] for r in rows] == ['0', '10']
    assert float(rows[0][3]) == round(want0, 5)
    # default: the control row without reads stays, as NA; the other values do not change
    (header, rows_na), = score.score([(t, c)], survival_fraction=0.5)
    assert rows_na[:2] == rows and rows_na[2] == ['chr1', '20', '30', 'NA']


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
    cfg.write_text('')
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
    cfg.write_text('')
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
    cfg.write_text('')
    seen = []
    monkeypatch.setattr(cli.subprocess, 'call', lambda cmd, **k: seen.append(cmd) or 0)
    monkeypatch.setattr(cli.subprocess, 'run', lambda cmd, **k: seen.append(cmd) or type('R', (), {'stdout': '1'})())
    monkeypatch.delenv('CONDENSE_CONTROLLER', raising=False)
    cli.main(['run', str(cfg), '-o', str(tmp_path / 'out'), '--slurm'])
    assert seen[0][0] == 'sbatch' and '--export=ALL,CONDENSE_CONTROLLER=1' in seen[0]


def _local_cmd(tmp_path, monkeypatch, extra, env):
    from condense_calc import cli
    cfg = tmp_path / 'config.yaml'
    cfg.write_text('')
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
    monkeypatch.setattr(cli, '_check', lambda src: None)
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
    monkeypatch.setattr(cli, '_check', lambda src: None)
    monkeypatch.setattr(reference, 'missing', lambda g, r=None: 'download hg38 into X')
    ids = iter(['111', '222'])
    seen = []
    monkeypatch.setattr(cli.subprocess, 'run',
                        lambda cmd, **k: seen.append(cmd) or type('R', (), {'stdout': next(ids)})())
    cli.main(['run', str(sh), '-o', str(tmp_path / 'out'), '--slurm'])
    ref, ctl = seen
    assert 'reference' in ref[-1] and ref[ref.index('-c') + 1] == '8'
    assert '--dependency=afterok:111' in ctl and ctl[ctl.index('-c') + 1] == '1'


def test_settings_resolve_and_validate():
    from condense_calc import settings
    S = settings.resolve({})
    assert S == {**settings.DEFAULTS, 'park': False}
    S = settings.resolve({'reads': {'dedup': 'true'}, 'align_threads': 24, 'samples': 'x.tsv', 'qc': 'OFF'})
    assert S['reads'] == {**settings.DEFAULTS['reads'], 'dedup': True}      # merged one level deep
    assert S['threads'] == 24 and S['qc'] == 'off'                          # earlier name of threads
    assert settings.resolve({'park': 'true'})['park'] is True                # as Snakemake's --config park=true
    with pytest.raises(settings.SettingsError, match="did you mean 'reads'"):
        settings.resolve({'read': {'min_mapq': 10}})
    with pytest.raises(settings.SettingsError, match="reads.min_mapqq"):
        settings.resolve({'reads': {'min_mapqq': 10}})
    with pytest.raises(settings.SettingsError, match='qc_thresholds.min_replicate_rh'):
        settings.resolve({'qc_thresholds': {'min_replicate_rh': 0.7}})
    with pytest.raises(settings.SettingsError, match='qc must be one of'):
        settings.resolve({'qc': 'full'})


def test_stage_defaults_from_settings():
    import argparse
    from condense_calc import cli, settings
    a = argparse.Namespace(park=False, mm=None, min_mapq=None, dedup=None, missing_as_nan=True)
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


def test_park_reaches_the_workflow(tmp_path, monkeypatch):
    # a config run gets --config park=true; a sheet run has park: true in condense.yaml
    cmd = _local_cmd(tmp_path, monkeypatch, ['--cores', '1', '--park'], {})
    assert cmd[cmd.index('--config') + 1:cmd.index('--config') + 2] == ['park=true']


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


def test_park_preset_and_explicit_keys():
    import argparse
    from condense_calc import cli, settings
    assert settings.resolve({})['reads']['min_mapq'] == 10
    park = settings.resolve({'park': True})
    assert park['reads'] == {'mm': 10, 'min_mapq': 0, 'dedup': False}
    assert settings.resolve({'park': True, 'reads': {'min_mapq': 5}})['reads']['min_mapq'] == 5    # explicit wins
    a = argparse.Namespace(park=True, min_mapq=None, mm=None)
    cli._stage_defaults(a)
    assert (a.min_mapq, a.mm) == (0, 10)

"""Settings of a run: the defaults, the Park preset, and their validation.

One source of truth for the workflow (Snakefile), `condense run` and the stage subcommands.

`park: true` (`condense run --park`, `condense-calc <stage> --park`) reproduces the processing that made
the GEO tables of Park et al. 2025 (GSE252941) with Park's scripts. PARK holds the settings that differ
from the defaults; in Park mode the stages also switch to Park's own behaviour where that has no setting.
Keys the user gives explicitly override both.
"""
import difflib
import re

from .reference import VERSION

DEFAULTS = {
    'park': False,
    'trim_to': 50,                      # reads are trimmed to this length before alignment
    'threads': 8,                       # bowtie2 threads per library
    'bin_sizes': [1000, 10000],
    'reads': {'mm': 10, 'min_mapq': 10, 'dedup': False},
    'score': {'missing_as_nan': True},  # bins without input reads: NA rows, so every table has the full bin grid
    'mscale': {'Ca': 5},                # molecule-number scale per agent (others 1)
    'fit': {'bin_size': 10000, 'min_tnum': 5, 'model': 'sigmoid', 'method': 'evolution', 'min_rsq': 0.5},
    'bowtie2': 'bowtie2',               # 'bowtie2' (PATH: this environment's), a version (installed once), or a path
    'fastq_check': 'lite',              # before alignment: stop on SRA Lite reads (fastq.py)
    'bdata': False,                     # deep series: also sliding-window areas (bindata_ver3.py); no table uses them
    'qc': 'default',
    'qc_thresholds': {},
    'max_runtime': 3 * 24 * 60,         # minutes; the longest time limit a retried step may ask for
}
PARK = {
    'reads': {'min_mapq': 0},         # Park's MAPQ filter is commented out (bincount_ver3.py:111-114)
    'score': {'missing_as_nan': False},   # NCPscore_ver4.py drops bins without input reads
    'bowtie2': '2.4.1',                 # the version in Park et al.'s Reporting Summary
    'fastq_check': 'park',              # also stop on read names not in the form Park aligned
}

# Park mode: the reference Park et al. used in place of a preset (reference.py; mm10 is already theirs)
PARK_GENOMES = {'hg38': 'hg38_ucsc'}

# Park mode, when `chroms` is not given: Park's driver (batch_submit_ver7.py:184-189) scores the female cells
# GM (GM12878) and mCD8T (WT, DFMO, ODCKO) without chrY, the others (H1, E14) with it
PARK_NO_CHRY = ('GM', 'mCD8T')


def park_chroms(label, chroms):
    """chroms without chrY if label is GM or mCD8T: a cell name, or a series name that starts with one
    (GM_NCP_sp, mCD8T_WT_NCP_sp, mCD8T:DFMO)."""
    if re.split(r'[_:-]', label)[0] in PARK_NO_CHRY:
        return [c for c in chroms if c not in ('chrY', 'Y')]
    return list(chroms)


def cell_chroms(label, given, default, park=False):
    """Chromosomes of a cell (config) or series (sample sheet): `chroms` if given, which always wins; else
    default(), the genome's; in Park mode without chrY for GM and mCD8T (park_chroms). Returns (chroms, warning):
    in Park mode a warning if a cell keeps chrY by default although its name suggests GM or mCD8T, because
    park_chroms goes by the first part of the name only (GM12878_NCP_sp keeps chrY)."""
    if given:
        return list(given), None
    chroms = park_chroms(label, default()) if park else list(default())
    if park and re.search(r'GM|CD8', label) and any(c in ('chrY', 'Y') for c in chroms):
        return chroms, (f"Park mode: {label} keeps chrY. Park et al. scored GM and mCD8T without it, and Park mode "
                        f"knows them by the first part of the name (GM_..., mCD8T_...); if {label} is one of them, "
                        "give chroms without chrY")
    return chroms, None


# keys that say what to run, not how: written by `condense` or given in a config file
INPUTS = ('sheet', 'genome', 'genomes', 'nucleosomes', 'samples', 'titration_dir', 'cells', 'chroms', 'refs',
          'align_threads')
BOOLS = ('park', 'bdata', 'reads.dedup', 'score.missing_as_nan')
INTS = ('trim_to', 'threads', 'max_runtime', 'reads.mm', 'reads.min_mapq', 'fit.bin_size', 'fit.min_tnum')
FLOATS = ('fit.min_rsq',)               # and the values of the NAMED dicts
CHOICES = {'qc': ('off', 'default', 'verbose'), 'fastq_check': ('off', 'lite', 'park'),
           'fit.model': ('sigmoid', 'hill'), 'fit.method': ('evolution', 'curve_fit')}
NAMED = ('mscale', 'qc_thresholds')     # dicts whose keys are names (agents, QC thresholds)
CELL_KEYS = ('genome', 'chroms')        # of each entry under `cells:` in a config


class SettingsError(ValueError):
    pass


def _unknown(key, known):
    near = difflib.get_close_matches(key, known, n=1)
    hint = f"; did you mean '{near[0]}'?" if near else f'; known: {", ".join(known)}'
    return SettingsError(f"unknown setting '{key}'{hint}")


def _bool(v, key):
    if isinstance(v, bool):
        return v
    if str(v).lower() in ('true', 'yes', '1'):
        return True
    if str(v).lower() in ('false', 'no', '0'):
        return False
    raise SettingsError(f'{key} must be true or false, not {v!r}')


def _number(v, key, kind=float):
    """v if it is an int (kind int) or a number (kind float); numbers written as text, as Snakemake's --config
    gives nested values, are converted."""
    if not isinstance(v, bool) and (isinstance(v, int) or (kind is float and isinstance(v, float))):
        return v
    if isinstance(v, str):
        try:
            return kind(v)
        except ValueError:
            pass
    raise SettingsError(f'{key} must be {"an integer" if kind is int else "a number"}, not {v!r}')


def _slot(out, key):
    """(dict, key) that hold a setting, also a dotted one such as 'reads.min_mapq'."""
    sec, _, k = key.rpartition('.')
    return (out[sec] if sec else out), k


def is_park(config):
    return _bool((config or {}).get('park', False), 'park')


def resolve(config=None, park=None):
    """Effective settings: DEFAULTS, then PARK in Park mode, then the user's keys (dicts merged one level deep).

    config: the user's settings, or a whole workflow config (its INPUTS keys are skipped). park, if given,
    overrides config['park']. Unknown keys and invalid values raise SettingsError."""
    user = dict(config or {})
    if park is not None:
        user['park'] = park
    if 'align_threads' in user and 'threads' not in user:      # the earlier name of `threads`
        user['threads'] = user['align_threads']
    for key in user:
        if key not in DEFAULTS and key not in INPUTS:
            raise _unknown(key, sorted(DEFAULTS))
    for c, v in (user.get('cells') or {}).items():
        for k in (v if isinstance(v, dict) else ()):
            if k not in CELL_KEYS:
                raise _unknown(f'cells.{c}.{k}', [f'cells.{c}.{x}' for x in CELL_KEYS])
    mode = is_park(user)
    out = {}
    for key, default in DEFAULTS.items():
        preset = PARK.get(key) if mode else None
        given = user.get(key)
        if isinstance(default, dict):
            if given is not None and not isinstance(given, dict):
                raise SettingsError(f'{key} must be a mapping, e.g. {key}: {default}')
            given = given or {}
            if key not in NAMED:
                for k in given:
                    if k not in default:
                        raise _unknown(f'{key}.{k}', [f'{key}.{d}' for d in default])
            out[key] = {**default, **(preset or {}), **given}
        else:
            out[key] = given if key in user else (preset if preset is not None else default)
    out['park'] = mode
    out['bowtie2'] = str(out['bowtie2'])        # YAML reads 2.5 as a number
    if re.fullmatch(r'[\d.]+', out['bowtie2']) and not VERSION.fullmatch(out['bowtie2']):
        raise SettingsError(f"bowtie2 must be a full version such as 2.4.1, not {out['bowtie2']} "
                            "(YAML reads 2.10 as the number 2.1: write '2.10.0' or quote it)")
    for key in BOOLS:
        d, k = _slot(out, key)
        d[k] = _bool(d[k], key)
    for key, choices in CHOICES.items():
        d, k = _slot(out, key)
        d[k] = 'off' if d[k] is False else str(d[k]).lower()     # YAML reads a bare off as false
        if d[k] not in choices:
            raise SettingsError(f'{key} must be one of {", ".join(choices)}, not {d[k]!r}')
    for key in INTS + FLOATS:
        d, k = _slot(out, key)
        d[k] = _number(d[k], key, int if key in INTS else float)
    for sec in NAMED:
        out[sec] = {k: _number(v, f'{sec}.{k}') for k, v in out[sec].items()}
    b = out['bin_sizes']
    if not isinstance(b, list) or not b:
        raise SettingsError(f'bin_sizes must be a list of bin sizes in bp, e.g. [1000, 10000], not {b!r}')
    out['bin_sizes'] = [_number(x, 'bin_sizes', int) for x in b]
    if min(out['bin_sizes']) <= 0:
        raise SettingsError(f'bin_sizes must be positive, not {b!r}')
    from .qc import THRESHOLDS
    for k in out['qc_thresholds']:
        if k not in THRESHOLDS:
            raise _unknown(f'qc_thresholds.{k}', [f'qc_thresholds.{t}' for t in THRESHOLDS])
    return out


def describe(section, key=None):
    """'default X' or 'default X; Park mode Y', for help texts."""
    d = DEFAULTS[section] if key is None else DEFAULTS[section][key]
    p = PARK.get(section) if key is None else PARK.get(section, {}).get(key)
    return f'default {d}' + (f'; --park: {p}' if p is not None and p != d else '')

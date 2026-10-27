"""Settings of a run: the defaults, the Park preset, and their validation.

One source of truth for the workflow (Snakefile), `condense run` and the stage subcommands.

`park: true` (`condense run --park`, `condense-calc <stage> --park`) reproduces the processing that made
the GEO tables of Park et al. 2025 (GSE252941) with Park's scripts. PARK holds the settings that differ
from the defaults; in Park mode the stages also switch to Park's own behaviour where that has no setting.
Keys the user gives explicitly override both.
"""
import difflib
import re

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


# keys that say what to run, not how: written by `condense` or given in a config file
INPUTS = ('sheet', 'genome', 'genomes', 'nucleosomes', 'samples', 'titration_dir', 'cells', 'chroms', 'refs',
          'align_threads')
BOOLS = ('park', 'bdata', 'reads.dedup', 'score.missing_as_nan')
CHOICES = {'qc': ('off', 'default', 'verbose'), 'fastq_check': ('off', 'lite', 'park')}
NAMED = ('mscale', 'qc_thresholds')     # dicts whose keys are names (agents, QC thresholds)


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
    for key in BOOLS:
        sec, _, k = key.rpartition('.')
        d = out[sec] if sec else out
        d[k] = _bool(d[k], key)
    for key, choices in CHOICES.items():
        out[key] = 'off' if out[key] is False else str(out[key]).lower()     # YAML reads a bare off as false
        if out[key] not in choices:
            raise SettingsError(f'{key} must be one of {", ".join(choices)}, not {out[key]!r}')
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

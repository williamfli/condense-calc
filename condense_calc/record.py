"""The run record, run_info.json in the output folder: what made the tables.

Mode, effective settings, condense-calc version and commit, bowtie2 version, references (with their
fingerprints, reference.py), chromosomes per cell, and per library the bowtie2 version (line and number) and the
FASTQ check (logs/align/). In Park mode, park_warnings() lists what keeps the tables from being Park's.
"""
import json
import shutil
import subprocess
from pathlib import Path

from . import __version__, reference

PKG_PARENT = Path(__file__).resolve().parents[1]


def code_version():
    """condense-calc's version, plus the commit of a git checkout ('-dirty' with uncommitted changes)."""
    git = ['git', '-C', str(PKG_PARENT)]
    try:
        commit = subprocess.run(git + ['rev-parse', '--short=12', 'HEAD'], capture_output=True, text=True,
                                check=True).stdout.strip()
        dirty = subprocess.run(git + ['status', '--porcelain', '--untracked-files=no'], capture_output=True,
                               text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return __version__
    return f'{__version__} ({commit}{"-dirty" if dirty else ""})'


def genome_info(fasta, index):
    fai = Path(f'{fasta}.fai')
    rows = [l for l in fai.read_text().splitlines() if l.strip()] if fai.exists() else []
    return {'fasta': str(fasta), 'bowtie2_index': str(index), 'sequences': len(rows),
            'fingerprint': reference.fingerprint(fai) if rows else None,
            'preset': reference.identify(fai) if rows else None}


def park_warnings(settings, genomes):
    """genomes: {name: {'fasta', 'bowtie2_index'}} used by the run."""
    if not settings['park']:
        return []
    out = []
    for g, v in genomes.items():
        preset = genome_info(v['fasta'], v['bowtie2_index'])['preset']
        if preset not in reference.PARK:
            out.append(f"Park mode, but genome {g} ({v['fasta']}) is not a reference Park et al. used "
                       f"({' or '.join(reference.PARK)}, by contig fingerprint): the tables will differ from GEO's")
    return out


def _first_line(path):
    return Path(path).read_text().splitlines()[0].strip() if Path(path).exists() else None


def write(path, settings, genomes, bt2_dir, chroms, names, logs='logs/align'):
    try:
        bt2_version = reference.bowtie2_version(bt2_dir)
    except (OSError, subprocess.CalledProcessError):
        bt2_version = None
    libs = {}
    for n in names:
        check, version = Path(logs) / f'{n}.fastq_check.json', Path(logs) / f'{n}.bowtie2_version'
        libs[n] = {'bowtie2': _first_line(version),
                   'bowtie2_version': reference.parse_version(version.read_text()) if version.exists() else None,
                   'fastq_check': json.loads(check.read_text()) if check.exists() else None}
    rec = {'mode': 'park' if settings['park'] else 'default',
           'condense_calc': code_version(),
           'settings': settings,
           'bowtie2': {'setting': settings['bowtie2'], 'version': bt2_version,
                       'folder': bt2_dir or str(Path(shutil.which('bowtie2') or '').parent)},
           'references': {g: genome_info(v['fasta'], v['bowtie2_index']) for g, v in genomes.items()},
           'chroms': chroms,
           'libraries': libs,
           'warnings': park_warnings(settings, genomes)}
    Path(path).write_text(json.dumps(rec, indent=1) + '\n')

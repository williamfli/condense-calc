"""Genome references: FASTA (+ .fai) and bowtie2 index, downloaded on first use (get_reference.sh).

Presets, each identified by a fingerprint of its sequence names and lengths (md5 of the first two
columns of the .fai), which is checked whenever a preset is used:
  hg38       NCBI GRCh38 no-alt analysis set: no alt haplotypes; chrY PARs and duplicate arrays masked
  mm10       ENCODE mm10 no-alt analysis set (ENCFF309GLL): 66 sequences, chrY PAR masked. Park et al.
             used the same index (4DN 4DNFI2493SDN), so Park mode keeps it
  hg38_ucsc  UCSC hg38.fa as distributed: 455 sequences with 261 alt haplotypes, nothing masked; Park
             et al.'s human reference, used for hg38 in Park mode (settings.PARK_GENOMES)
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

PKG = Path(__file__).resolve().parent
PRESETS = {
    'hg38': {'name': 'NCBI GRCh38 no-alt analysis set', 'contigs': 195,
             'fingerprint': '9b996e0d2eabf4f49ef31793c11d2f45'},
    'mm10': {'name': 'ENCODE mm10 no-alt analysis set (ENCFF309GLL)', 'contigs': 66,
             'fingerprint': '1a3951070d28838887ce0c77dc27c2e3'},
    'hg38_ucsc': {'name': 'UCSC hg38.fa with alt haplotypes', 'contigs': 455,
                  'fingerprint': '65b2bdfb1382e0a8a219be9c04a8f1a7'},
}
# sequence sets met before under a preset's name
OTHERS = {'a0b9649980a9063a8d14ff2bfd7c202c': "UCSC mm10 chr1-19, X, Y, M (22 sequences; condense-calc's "
                                              'mm10 until 2026-10)'}


CANONICAL = re.compile(r'(chr)?([0-9]+|X|Y)')
MITO = {'chrM', 'chrMT', 'M', 'MT', 'mitochondrion'}


def default_chroms(fai):
    """Chromosomes processed when the settings name none, from a FASTA index (.fai).

    hg38, mm10 and assemblies like them (chromosomes named 1..N, X, Y, with or without 'chr', making up
    >= 90% of the sequence): exactly those, in FASTA order, as in the GEO run. Any other genome
    (e.g. dm6 chr2L, yeast chrI-chrXVI, NCBI NC_ names): every sequence at least 1% as long as the
    longest one, except the mitochondrial genome; this leaves out small unplaced contigs.
    """
    rows = [l.split('\t')[:2] for l in Path(fai).read_text().splitlines() if l.strip()]
    if not rows:
        raise ValueError(f'{fai} lists no sequences')
    size = {n: int(length) for n, length in rows}
    canon = [n for n, _ in rows if CANONICAL.fullmatch(n)]
    if canon and sum(size[n] for n in canon) >= 0.9 * sum(size.values()):
        return canon
    longest = max(size.values())
    keep = [n for n, _ in rows if n not in MITO and size[n] >= 0.01 * longest]
    if not keep:
        raise ValueError(f'{fai}: no chromosomes besides the mitochondrial genome')
    return keep


def fingerprint(fai):
    """md5 of the sequence names and lengths of a FASTA index (= `cut -f1,2 FAI | md5sum`)."""
    rows = [l.split('\t')[:2] for l in Path(fai).read_text().splitlines() if l.strip()]
    return hashlib.md5(''.join(f'{n}\t{L}\n' for n, L in rows).encode()).hexdigest()


def identify(fai):
    """The preset a FASTA index belongs to, or None."""
    fp = fingerprint(fai)
    return next((g for g, p in PRESETS.items() if p['fingerprint'] == fp), None)


def _check_preset(genome, fasta):
    """Stop if a preset's folder holds another sequence set, e.g. an older install under the same name."""
    fai = Path(f'{fasta}.fai')
    if genome not in PRESETS or not fai.exists():
        return
    fp = fingerprint(fai)
    if fp != PRESETS[genome]['fingerprint']:
        n = sum(1 for l in fai.read_text().splitlines() if l.strip())
        held = OTHERS.get(fp, f'{n} sequences, fingerprint {fp}')
        raise SystemExit(f'{fasta.parent} holds another reference ({held}), not {genome} = '
                         f'{PRESETS[genome]["name"]}.\nMove it away (e.g. mv {fasta.parent} {fasta.parent}.old) '
                         f'and run again: {genome} is then downloaded. Or give that FASTA as --genome.')


def refs_dir(arg=None):
    return Path(arg or os.environ.get('CONDENSE_REFS') or Path.home() / '.cache' / 'condense-calc').expanduser()


VERSION = re.compile(r'\d+(\.\d+)+')


def bowtie2_dir(spec, refs=None):
    """Folder of the bowtie2 to run, or None for the one on PATH. spec (settings `bowtie2`): 'bowtie2' (PATH),
    a version (a conda env <refs>/envs/bowtie2-<version>, made by ensure_bowtie2), or a folder or file path."""
    spec = str(spec)
    if spec == 'bowtie2':
        return None
    if VERSION.fullmatch(spec):
        return str(refs_dir(refs) / 'envs' / f'bowtie2-{spec}' / 'bin')
    p = Path(spec).expanduser().resolve()
    if p.is_dir():
        return str(p / 'bin' if (p / 'bin' / 'bowtie2').exists() else p)
    return str(p.parent)



def missing_bowtie2(spec, refs=None):
    """What ensure_bowtie2() would still have to do: None if that bowtie2 is there, else a short description."""
    d = bowtie2_dir(spec, refs)
    if d is None or (Path(d) / 'bowtie2').exists():
        return None
    if VERSION.fullmatch(str(spec)):
        return f'install bowtie2 {spec} (a conda env in {Path(d).parent})'
    raise SystemExit(f'bowtie2: no bowtie2 in {d} (settings bowtie2: {spec})')


def bowtie2_version(folder=None):
    """Version of the bowtie2 in folder (None: on PATH). Its wrapper is a Perl script, so the folder goes first on
    PATH, as the align step does. Some conda builds print no version (bioconda 2.4.1: 'bowtie2-align-s version ');
    then the environment's conda-meta record says it."""
    env = dict(os.environ, PATH=f'{folder}{os.pathsep}{os.environ.get("PATH", "")}') if folder else None
    out = subprocess.run(['bowtie2', '--version'], env=env, capture_output=True, text=True, check=True).stdout
    m = re.match(r'\S+ version (\S+)', out)
    if m:
        return m.group(1)
    where = Path(folder) if folder else Path(shutil.which('bowtie2') or '.').parent
    meta = sorted(where.parent.glob('conda-meta/bowtie2-[0-9]*.json'))
    return json.loads(meta[-1].read_text())['version'] if meta else None


def ensure_bowtie2(spec, refs=None):
    """Folder of the bowtie2 to run (None: PATH). A version is installed once, with mamba or conda, as
    <refs>/envs/bowtie2-<version>; only bowtie2 runs from it (samtools stays this environment's)."""
    todo = missing_bowtie2(spec, refs)
    d = bowtie2_dir(spec, refs)
    if todo:
        conda = (os.environ.get('MAMBA_EXE') or os.environ.get('CONDA_EXE') or shutil.which('mamba')
                 or shutil.which('conda'))
        if not conda:
            raise SystemExit(f'bowtie2 {spec}: mamba or conda is needed to install it; or set bowtie2: to the '
                             'folder of an installed bowtie2')
        env = Path(d).parent
        print(f'Installing bowtie2 {spec} into {env}; this happens once.', flush=True)
        env.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([conda, 'create', '-y', '-q', '-p', str(env), '-c', 'conda-forge', '-c', 'bioconda',
                        f'bowtie2={spec}'], check=True)
    if d is not None and VERSION.fullmatch(str(spec)) and bowtie2_version(d) != str(spec):
        raise SystemExit(f'{d}/bowtie2 is version {bowtie2_version(d)}, not {spec}')
    return d


def _index_ok(prefix):
    return len(list(Path(prefix).parent.glob(Path(prefix).name + '.*.bt2'))) + \
        len(list(Path(prefix).parent.glob(Path(prefix).name + '.*.bt2l'))) >= 6


def locate(genome, refs=None):
    """(name, fasta, bowtie2 index prefix) of a genome name or FASTA path, without fetching anything."""
    if genome in PRESETS:
        d = refs_dir(refs) / genome
        return genome, d / f'{genome}.fa', d / genome
    fasta = Path(genome).expanduser().resolve()
    if not fasta.exists():
        raise SystemExit(f'--genome: {genome} is neither {", ".join(PRESETS)} nor a FASTA file')
    index = fasta.with_suffix('') if fasta.suffix in ('.fa', '.fasta', '.fna') else fasta
    return index.name, fasta, index


def missing(genome, refs=None):
    """What ensure() would still have to do: None if the genome is ready, else a short description."""
    name, fasta, index = locate(genome, refs)
    _check_preset(genome, fasta)
    if Path(f'{fasta}.fai').exists() and _index_ok(index):
        return None
    if genome in PRESETS:
        return f'download {genome} (FASTA and bowtie2 index) into {fasta.parent}'
    return f'build the bowtie2 index {index} (large genomes take hours)'


def ensure(genome, refs=None, threads=1):
    """Return (name, fasta, bowtie2 index prefix) for a genome name or a FASTA path."""
    if genome in PRESETS:
        d = refs_dir(refs) / genome
        fasta, index = d / f'{genome}.fa', d / genome
        _check_preset(genome, fasta)
        if not (Path(f'{fasta}.fai').exists() and _index_ok(index)):
            print(f'Downloading {genome} (FASTA and bowtie2 index) into {d}; this happens once.', flush=True)
            d.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(['bash', str(PKG / 'get_reference.sh'), str(d.parent), str(threads), genome], check=True)
            _check_preset(genome, fasta)
        return genome, str(fasta), str(index)
    fasta = Path(genome).expanduser().resolve()
    if not fasta.exists():
        raise SystemExit(f'--genome: {genome} is neither {", ".join(PRESETS)} nor a FASTA file')
    index = fasta.with_suffix('') if fasta.suffix in ('.fa', '.fasta', '.fna') else fasta
    if not Path(f'{fasta}.fai').exists():
        subprocess.run(['samtools', 'faidx', str(fasta)], check=True)
    if not _index_ok(index):
        print(f'Building the bowtie2 index {index} (large genomes take hours).', flush=True)
        subprocess.run(['bowtie2-build', '--threads', str(threads), str(fasta), str(index)], check=True)
    return index.name, str(fasta), str(index)

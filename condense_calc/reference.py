"""Genome references: FASTA (+ .fai) and bowtie2 index, downloaded on first use."""
import os
import re
import subprocess
from pathlib import Path

PKG = Path(__file__).resolve().parent
KNOWN = ('hg38', 'mm10')


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


def refs_dir(arg=None):
    return Path(arg or os.environ.get('CONDENSE_REFS') or Path.home() / '.cache' / 'condense-calc').expanduser()


def _index_ok(prefix):
    return len(list(Path(prefix).parent.glob(Path(prefix).name + '.*.bt2'))) + \
        len(list(Path(prefix).parent.glob(Path(prefix).name + '.*.bt2l'))) >= 6


def locate(genome, refs=None):
    """(name, fasta, bowtie2 index prefix) of a genome name or FASTA path, without fetching anything."""
    if genome in KNOWN:
        d = refs_dir(refs) / genome
        return genome, d / f'{genome}.fa', d / genome
    fasta = Path(genome).expanduser().resolve()
    if not fasta.exists():
        raise SystemExit(f'--genome: {genome} is neither {" nor ".join(KNOWN)} nor a FASTA file')
    index = fasta.with_suffix('') if fasta.suffix in ('.fa', '.fasta', '.fna') else fasta
    return index.name, fasta, index


def missing(genome, refs=None):
    """What ensure() would still have to do: None if the genome is ready, else a short description."""
    name, fasta, index = locate(genome, refs)
    if Path(f'{fasta}.fai').exists() and _index_ok(index):
        return None
    if genome in KNOWN:
        return f'download {genome} (FASTA and bowtie2 index) into {fasta.parent}'
    return f'build the bowtie2 index {index} (large genomes take hours)'


def ensure(genome, refs=None, threads=1):
    """Return (name, fasta, bowtie2 index prefix) for a genome name or a FASTA path."""
    if genome in KNOWN:
        d = refs_dir(refs) / genome
        fasta, index = d / f'{genome}.fa', d / genome
        if not (Path(f'{fasta}.fai').exists() and _index_ok(index)):
            print(f'Downloading {genome} (FASTA and bowtie2 index) into {d}; this happens once.', flush=True)
            d.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(['bash', str(PKG / 'get_reference.sh'), str(d.parent), str(threads), genome], check=True)
        return genome, str(fasta), str(index)
    fasta = Path(genome).expanduser().resolve()
    if not fasta.exists():
        raise SystemExit(f'--genome: {genome} is neither {" nor ".join(KNOWN)} nor a FASTA file')
    index = fasta.with_suffix('') if fasta.suffix in ('.fa', '.fasta', '.fna') else fasta
    if not Path(f'{fasta}.fai').exists():
        subprocess.run(['samtools', 'faidx', str(fasta)], check=True)
    if not _index_ok(index):
        print(f'Building the bowtie2 index {index} (large genomes take hours).', flush=True)
        subprocess.run(['bowtie2-build', '--threads', str(threads), str(fasta), str(index)], check=True)
    return index.name, str(fasta), str(index)

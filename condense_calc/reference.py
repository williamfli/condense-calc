"""Genome references: FASTA (+ .fai) and bowtie2 index, downloaded on first use."""
import os
import subprocess
from pathlib import Path

PKG = Path(__file__).resolve().parent
KNOWN = ('hg38', 'mm10')


def refs_dir(arg=None):
    return Path(arg or os.environ.get('CONDENSE_REFS') or Path.home() / '.cache' / 'condense-calc').expanduser()


def _index_ok(prefix):
    return len(list(Path(prefix).parent.glob(Path(prefix).name + '.*.bt2'))) + \
        len(list(Path(prefix).parent.glob(Path(prefix).name + '.*.bt2l'))) >= 6


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

"""`condense fetch SRR...`: the FASTQ files of SRA runs as Park et al. aligned them.

From NCBI's SRA Normalized data (real base qualities; ENA serves some GSE252941 runs as SRA Lite), with the
instrument's read names '@<spot name> <mate>:N:0:<index>' (bowtie2 places multi-mapping reads by name). These
pass `fastq_check: park`: `condense fetch SRR... -o reads/`, then `condense samples.tsv --park` with the
titration at full precision, reproduces Park's GEO tables.
Per run: prefetch (--verify yes; up to 4 tries, waiting 2, 4, 6 min), fasterq-dump, a check that both mates
hold the run's spot count and real qualities, gzip (pigz if installed), then the .sra and temporary files are
deleted. Writes OUT/<run>_1.fastq.gz and _2.fastq.gz; a run whose two files exist is skipped.
"""
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from . import fastq

DEFLINE = '@$sn $ri:N:0:$sg'
TRIES, WAIT = 4, 120


class FetchError(RuntimeError):
    pass


def _tool(name):
    path = shutil.which(name)
    if not path:
        raise FetchError(f'{name} not found: install sra-tools (conda install -c bioconda sra-tools)')
    return path


def _reads(path):
    with open(path, 'rb') as f:
        return sum(chunk.count(b'\n') for chunk in iter(lambda: f.read(1 << 24), b'')) // 4


def fetch(run, out, threads=1, sleep=time.sleep):
    """OUT/<run>_1.fastq.gz and OUT/<run>_2.fastq.gz; raises FetchError."""
    out = Path(out)
    final = [out / f'{run}_{m}.fastq.gz' for m in (1, 2)]
    if all(f.exists() for f in final):
        print(f'{run}: already there', flush=True)
        return final
    tmp = out / f'tmp_{run}'
    tmp.mkdir(parents=True, exist_ok=True)
    sra = tmp / run / f'{run}.sra'
    prefetch = [_tool('prefetch'), '--verify', 'yes', '--max-size', '200G', '-O', str(tmp), run]
    for t in range(1, TRIES + 1):
        if subprocess.run(prefetch).returncode == 0:
            break
        if t == TRIES:
            raise FetchError(f'{run}: prefetch failed {TRIES} times')
        print(f'{run}: prefetch try {t} of {TRIES} failed; next try in {WAIT * t // 60} min', file=sys.stderr, flush=True)
        sleep(WAIT * t)
    if not sra.exists():
        raise FetchError(f'{run}: prefetch gave no SRA Normalized {sra.name} (only SRA Lite?): '
                         f'{", ".join(p.name for p in sra.parent.glob("*")) or "nothing"}')
    p = subprocess.run([_tool('fasterq-dump'), '--split-files', '-e', str(threads), '-f', '-t', str(tmp / 'fqd'),
                        '-O', str(tmp), '--seq-defline', DEFLINE, '--qual-defline', '+', '-o', f'{run}.fastq',
                        str(sra)], capture_output=True, text=True)
    print(p.stdout + p.stderr, end='', file=sys.stderr, flush=True)
    if p.returncode:
        raise FetchError(f'{run}: fasterq-dump failed (exit {p.returncode})')
    m = re.search(r'spots read\s*:\s*([\d,]+)', p.stdout + p.stderr)
    spots = int(m.group(1).replace(',', '')) if m else None
    fqs = [tmp / f'{run}_{m}.fastq' for m in (1, 2)]
    reads = [_reads(f) if f.exists() else None for f in fqs]
    if spots is None or reads != [spots, spots]:
        raise FetchError(f'{run}: {reads} reads in the two files, but {spots} spots in the run')
    probs = fastq.problems([fastq.check_file(f, mate, 'park') for mate, f in enumerate(fqs, 1)])
    if probs:
        raise FetchError('\n'.join(probs))
    gz = [_tool('pigz'), '-p', str(threads)] if shutil.which('pigz') else ['gzip']
    for f, dst in zip(fqs, final):
        Path(f'{f}.gz').unlink(missing_ok=True)      # left by a job stopped while compressing; gzip would refuse
        subprocess.run(gz + [str(f)], check=True)
        Path(f'{f}.gz').rename(dst)
    shutil.rmtree(tmp)
    print(f'{run}: {spots:,} read pairs: {final[0]}, {final[1].name}', flush=True)
    return final

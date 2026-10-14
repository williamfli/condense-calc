"""Make a synthetic condense-seq experiment with a known answer, for testing the pipeline.

    python tests/simulate.py OUTDIR [--seed 0]

Writes into OUTDIR:
  ref/sim.fa (+ .fai, bowtie2 index ref/sim)   3 chromosomes, 600 kb, with a 2 kb segment
                                               copied into chr1 and chr2 (multi-mappers) and
                                               a 5 kb N gap in chr2 (bins without reads)
  fastq/<name>_R1.fastq.gz, _R2.fastq.gz       paired 50-bp reads, named as the GEO run named BAMs
  titration/SIM_NCP_sp.tsv                     tnum, concentration, true survival fraction
  samples.tsv                                  sample sheet for the workflow
  truth/nucleosomes.tsv                        every nucleosome: dyad, true c-half, rate

Model: nucleosomes every ~190 bp; nucleosome i stays soluble at agent concentration c with
p_i(c) = 1 / (1 + exp(rate * (c - chalf_i))), where chalf_i varies in 60-kb domains. A
supernatant library samples nucleosomes in proportion to p_i(c); every library is sequenced
to a fixed depth. The true score of nucleosome i at c is then -ln p_i(c).
"""
import argparse
import gzip
import subprocess
from pathlib import Path

import numpy as np

CHROMS = {'chr1': 300_000, 'chr2': 200_000, 'chrX': 100_000}
CONC = [0.0, 0.2, 0.35, 0.5, 0.65, 0.8]          # tnum 0..5 (mM)
RATE = 12.0                                     # per mM
SHALLOW_PAIRS, DEEP_PAIRS = 40_000, 150_000
DEEP_TNUMS = [0, 2, 4]
READ_LEN = 50
COMP = str.maketrans('ACGTN', 'TGCAN')


def genome(rng):
    seqs = {}
    for c, n in CHROMS.items():
        seqs[c] = rng.choice(list('ACGT'), size=n, p=[0.29, 0.21, 0.21, 0.29])
    rep = seqs['chr1'][50_000:52_000].copy()          # duplicated segment -> multi-mappers
    seqs['chr2'][120_000:122_000] = rep
    seqs['chr2'][150_000:155_000] = 'N'               # assembly gap
    return {c: ''.join(s) for c, s in seqs.items()}


def nucleosomes(rng, seqs):
    rows = []
    for c, s in seqs.items():
        pos = 100
        while pos < len(s) - 300:
            if 'N' not in s[pos - 100:pos + 100]:
                dom = 0.45 + 0.15 * np.sin(2 * np.pi * pos / 60_000) + rng.normal(0, 0.04)
                rows.append((c, pos, dom))
            pos += int(rng.integers(160, 221))
    return rows


def fragments(rng, nucs, weights, n_pairs, seqs):
    """Sample fragments (chrom, start, length) from nucleosomes with the given weights."""
    p = np.asarray(weights) / np.sum(weights)
    pick = rng.choice(len(nucs), size=n_pairs, p=p)
    kind = rng.random(n_pairs)
    out = []
    for i, k in zip(pick, kind):
        c, dyad, _ = nucs[i]
        if k < 0.85:                                   # mono-nucleosome
            L = int(np.clip(rng.normal(147, 6), 120, 200))
        elif k < 0.95:                                 # sub-nucleosomal
            L = int(rng.integers(90, 121))
        else:                                          # di-nucleosome
            L = int(rng.integers(300, 360))
        centre = dyad + int(round(rng.normal(0, 3)))
        start = centre - L // 2
        if start < 0 or start + L > len(seqs[c]):
            continue
        out.append((c, start, L))
    return out


def mutate(rng, seq, rate=0.005):
    s = list(seq)
    for j in np.flatnonzero(rng.random(len(s)) < rate):
        s[j] = rng.choice([b for b in 'ACGT' if b != s[j]])
    return ''.join(s)


def write_fastq(path_base, frags, seqs, rng):
    q = 'I' * READ_LEN
    with gzip.open(f'{path_base}_R1.fastq.gz', 'wt') as f1, gzip.open(f'{path_base}_R2.fastq.gz', 'wt') as f2:
        for n, (c, start, L) in enumerate(frags):
            frag = seqs[c][start:start + L]
            left, right = frag[:READ_LEN], frag[-READ_LEN:].translate(COMP)[::-1]
            r1, r2 = (left, right) if rng.random() < 0.5 else (right, left)
            name = f'SIM:1:FC:1:{n // 1000}:{n % 1000}:{start}'
            f1.write(f'@{name} 1:N:0\n{mutate(rng, r1)}\n+\n{q}\n')
            f2.write(f'@{name} 2:N:0\n{mutate(rng, r2)}\n+\n{q}\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir')
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    out = Path(a.outdir)
    for d in ('ref', 'fastq', 'titration', 'truth'):
        (out / d).mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    seqs = genome(rng)
    with open(out / 'ref/sim.fa', 'w') as f:
        for c, s in seqs.items():
            f.write(f'>{c}\n')
            f.write('\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
    subprocess.run(['samtools', 'faidx', out / 'ref/sim.fa'], check=True)
    subprocess.run(['bowtie2-build', '-q', out / 'ref/sim.fa', out / 'ref/sim'], check=True)

    nucs = nucleosomes(rng, seqs)
    chalf = np.array([n[2] for n in nucs])
    with open(out / 'truth/nucleosomes.tsv', 'w') as f:
        f.write('chrom\tdyad\tchalf\trate\n')
        for c, d, h in nucs:
            f.write(f'{c}\t{d}\t{h:.5f}\t{RATE}\n')

    surv = []
    samples = ['cell\tsample\tagent\ttnum\trep\tdepth\tfastq1\tfastq2']
    for t, conc in enumerate(CONC):
        p = 1 / (1 + np.exp(RATE * (conc - chalf))) if conc > 0 else np.ones(len(nucs))
        surv.append(p.mean())
        for depth, n_pairs in (('shallow', SHALLOW_PAIRS), ('deep', DEEP_PAIRS)):
            if depth == 'deep' and t not in DEEP_TNUMS:
                continue
            name = f'SIM_NCP_sp_{t}_1rep' + ('_deep' if depth == 'deep' else '')
            frags = fragments(rng, nucs, p, n_pairs, seqs)
            write_fastq(out / 'fastq' / name, frags, seqs, rng)
            samples.append('\t'.join(['SIM', 'NCP', 'sp', str(t), '1', depth,
                                      f'fastq/{name}_R1.fastq.gz', f'fastq/{name}_R2.fastq.gz']))
    with open(out / 'samples.tsv', 'w') as f:
        f.write('\n'.join(samples) + '\n')
    with open(out / 'titration/SIM_NCP_sp.tsv', 'w') as f:
        f.write('tnum\tconc\tsurvival_fraction\n')
        for t, (c, s) in enumerate(zip(CONC, surv)):
            f.write(f'{t}\t{c}\t{s:.6f}\n')
    print(f'{len(nucs)} nucleosomes; survival fractions', [round(s, 3) for s in surv])


if __name__ == '__main__':
    main()

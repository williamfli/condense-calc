"""Fragments and low-MAPQ (< 10) fragments per 10 kb bin, to see which bins hold multi-mappers.

    python geo_test/mapq_bins.py BAM OUT.tsv

Counts each fragment once (forward read of a proper pair) at its centre, like bincount.
"""
import sys

import numpy as np
import pysam

bam, out = sys.argv[1], sys.argv[2]
f = pysam.AlignmentFile(bam); rows = []
for chrom in [f'chr{i}' for i in list(range(1,23))+['X','Y']]:
    n = f.get_reference_length(chrom)//10000 + 1
    tot = np.zeros(n, int); low = np.zeros(n, int)
    for r in f.fetch(chrom):
        if r.is_reverse or not r.is_proper_pair or r.is_unmapped or r.mate_is_unmapped: continue
        b = (r.reference_start + abs(r.template_length)//2)//10000
        if b >= n: continue
        tot[b] += 1; low[b] += r.mapping_quality < 10
    for i in range(n): rows.append(f'{chrom}\t{i*10000}\t{tot[i]}\t{low[i]}')
open(out,'w').write('Chromosome\tStart\ttot\tlowmapq\n'+'\n'.join(rows)+'\n')

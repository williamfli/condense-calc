"""Step 0, read level: how many reads pass each of Park's filters, and what the kept fragments look like.

    python scripts/01_bam_stats.py BAM_DIR RESULTS_DIR BIG_DIR

For each BAM (input = ID 0, supernatant = ID 6) it walks every record once and applies the
filters of condseq.fragments.read_bam in the same order (Park's read loop):
  forward-strand read (0x10 unset) -> mapped (0x4 unset) -> mate mapped (0x8 unset)
  -> proper pair (0x2 set) -> NM <= 10 -> on chr1-22, X, Y
The fragments left are exactly what `bincount` counts. For them it records MAPQ, fragment length
(|TLEN|) and, per 10 kb bin (the fragment centre as bincount places it), the total and the number
with MAPQ < 10 and with MAPQ <= 1.

Writes to RESULTS_DIR (small, kept in git):
  01_filter_cascade.tsv   reads left after each filter, per sample
  01_mapq_hist.tsv        kept fragments by MAPQ
  01_tlen_hist.tsv        kept fragments by |TLEN| (1-1000 bp; longer ones in the last row)
and to BIG_DIR: 01_bins_10kb_mapq.tsv.gz (one row per 10 kb bin and sample).
Env: condseq. ~10 min, 1 core, < 4 GB.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))   # condense-seq-pipeline/
from condseq.fragments import Fragments, _ref_span                # noqa: E402

CHROMS = [f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']]
SAMPLES = {'input (ID 0)': 'H1_NCP_PEG_0_2rep', 'supernatant (ID 6)': 'H1_NCP_PEG_6_2rep'}
BIN, MM = 10000, 10
STEPS = ['records', 'forward strand', 'mapped', 'mate mapped', 'proper pair', f'NM <= {MM}',
         'chr1-22, X, Y']


def stats(bam_path):
    bam = pysam.AlignmentFile(bam_path)
    n = dict.fromkeys(STEPS, 0)
    pos, end, tlen, mapq, chrom = [], [], [], [], []
    for r in bam.fetch(until_eof=True):
        n['records'] += 1
        f = r.flag
        if f & 0x10:
            continue
        n['forward strand'] += 1
        if f & 0x4:
            continue
        n['mapped'] += 1
        if f & 0x8:
            continue
        n['mate mapped'] += 1
        if not f & 0x2:
            continue
        n['proper pair'] += 1
        if r.has_tag('NM') and r.get_tag('NM') > MM:
            continue
        n[f'NM <= {MM}'] += 1
        if r.reference_name not in CHROMS:
            continue
        n['chr1-22, X, Y'] += 1
        pos.append(r.reference_start)
        end.append(r.reference_start + _ref_span(r.cigartuples))
        tlen.append(r.template_length)
        mapq.append(r.mapping_quality)
        chrom.append(r.reference_id)
    pos, end, tlen = (np.array(x, dtype=np.int64) for x in (pos, end, tlen))
    mapq, chrom = np.array(mapq), np.array(chrom)
    dyad = Fragments(pos, end, tlen).dyad_bincount()
    sizes = dict(zip(bam.references, bam.lengths))
    rows = []
    for name in CHROMS:
        sel = chrom == bam.get_tid(name)
        nb = sizes[name] // BIN + 1
        b = dyad[sel] // BIN
        b[b < 0] += nb                       # negative indices wrap, as in bincount
        q = mapq[sel]
        tot = np.bincount(b, minlength=nb)[:nb]
        low = np.bincount(b[q < 10], minlength=nb)[:nb]
        q01 = np.bincount(b[q <= 1], minlength=nb)[:nb]
        rows.append(pd.DataFrame({'Chromosome': name, 'Start': np.arange(nb) * BIN,
                                  'End': np.minimum(np.arange(1, nb + 1) * BIN, sizes[name]),
                                  'fragments': tot, 'mapq_lt10': low, 'mapq_le1': q01}))
    return n, mapq, np.abs(tlen), pd.concat(rows)


def main():
    bam_dir, res, big = (Path(a) for a in sys.argv[1:4])
    res.mkdir(parents=True, exist_ok=True)
    big.mkdir(parents=True, exist_ok=True)
    casc, mq, tl, bins = {}, {}, {}, []
    for label, name in SAMPLES.items():
        n, mapq, alen, b = stats(bam_dir / f'{name}.bam')
        casc[label] = n
        mq[label] = np.bincount(mapq, minlength=43)
        tl[label] = np.bincount(np.minimum(alen, 1001), minlength=1002)[1:]
        bins.append(b.assign(sample=name))
        print(name, n, flush=True)
    pd.DataFrame(casc).rename_axis('step').to_csv(res / '01_filter_cascade.tsv', sep='\t')
    pd.DataFrame(mq).rename_axis('mapq').to_csv(res / '01_mapq_hist.tsv', sep='\t')
    pd.DataFrame(tl, index=pd.RangeIndex(1, 1002, name='tlen')).to_csv(res / '01_tlen_hist.tsv', sep='\t')
    pd.concat(bins).to_csv(big / '01_bins_10kb_mapq.tsv.gz', sep='\t', index=False)


if __name__ == '__main__':
    main()

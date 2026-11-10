"""Chromatin-state annotations and cytogenetic bands for the Extended Data Fig. 2 QC panels (ed2.py).

`states` in the settings (`condense run --states`):
  none (default)  panels g-i are not drawn
  auto            a built-in annotation for the cell and genome (STATES below), downloaded once into
                  <refs>/annotations/<name>/ and checked against its md5; with it, the genome's cytogenetic bands
                  (CYTOBANDS) mask centromeres and heterochromatin in panel d
  FILE            a BED file (plain or gzipped): chrom, start, end, state name (more columns are ignored)
  {cell: ...}     per cell, each value as above
`cytoband: auto | FILE` gives the bands without states. Without bands, panel d masks assembly gaps (runs of N).

A cell matches a built-in annotation by the start of its name (GM, GM12878, E14, mESC; then the end of the name or
one of _ : - .), as the cells of a config or the series groups of a sample sheet are named. Other cells (H1, for
example, whose model in Park et al. is not published) need a BED file.
"""
import gzip
import re
from pathlib import Path

# ordered for the figures as in Park et al. (active -> repressed -> quiescent); (state name in the file, label)
STATES = {
    'GM12878_hg38': dict(
        cell=r'GM(12878)?', genome='hg38', file='GM12878_hg38/ENCFF140VIG.bed.gz',
        url='https://www.encodeproject.org/files/ENCFF140VIG/@@download/ENCFF140VIG.bed.gz',
        md5='30b7c8052eddf8c9fa1b36cb19bfdd55',
        source='ENCODE ENCFF140VIG (ENCSR310DFZ): ChromHMM 15-state model of GM12878, ENCODE4',
        order=[('Tss', 'Active TSS'), ('TssFlnk', 'TSS flank'), ('TssFlnkU', 'TSS flank up'),
               ('TssFlnkD', 'TSS flank down'), ('Biv', 'Bivalent'), ('EnhG1', 'Genic enhancer 1'),
               ('EnhG2', 'Genic enhancer 2'), ('Enh1', 'Enhancer 1'), ('Enh2', 'Enhancer 2'), ('Tx', 'Tx'),
               ('TxWk', 'Weak Tx'), ('ZNF/Rpts', 'ZNF / repeats'), ('ReprPC', 'Polycomb repressed'),
               ('Het', 'Heterochromatin'), ('Quies', 'Quiescent')]),
    'mESC_E14_mm10': dict(
        cell=r'E14|mESC', genome='mm10', file='mESC_E14_mm10/mESC_E14_12_segments.bed.gz',
        url=('https://raw.githubusercontent.com/guifengwei/ChromHMM_mESC_mm10/'
             '96506cbb5cd6be8e2bab3dc39300404f78cd8234/mESC_E14_12_segments.bed.gz'),
        md5='058c115e5d3b28956852664c50bd0814',
        source='Pintacuda, Wei et al., Mol Cell 2017: ChromHMM 12-state model of E14 mESC '
               '(github.com/guifengwei/ChromHMM_mESC_mm10, commit 96506cb)',
        order=[('7_ActivePromoter', 'Active promoter'), ('6_BivalentChromatin', 'Bivalent chromatin'),
               ('8_StrongEnhancer', 'Strong enhancer'), ('4_Enhancer', 'Enhancer'), ('11_WeakEnhancer', 'Weak enhancer'),
               ('9_TranscriptionTransition', 'Tx transition'), ('10_TranscriptionElongation', 'Tx elongation'),
               ('1_Insulator', 'Insulator'), ('5_RepressedChromatin', 'Repressed chromatin'),
               ('3_Heterochromatin', 'Heterochromatin'), ('2_Intergenic', 'Intergenic'),
               ('12_LowSignal/RepetitiveElements', 'Low signal / repeats')]),
}
CYTOBANDS = {
    'hg38': dict(file='cytoband_hg38/cytoBand.txt.gz', md5='6677dd50caf19b0c3626b99ce67ca073',
                 url='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz',
                 source='UCSC hg38 cytoBand'),
}
GENOME_BUILDS = {'hg38': 'hg38', 'hg38_ucsc': 'hg38', 'grch38': 'hg38', 'mm10': 'mm10', 'grcm38': 'mm10'}


def build_of(genome, fasta=None):
    """'hg38', 'mm10' or None: from the genome's name, else from its FASTA index (reference presets)."""
    b = GENOME_BUILDS.get(str(genome).lower())
    if b is None and fasta and Path(f'{fasta}.fai').exists():
        from . import reference
        b = GENOME_BUILDS.get(reference.identify(f'{fasta}.fai') or '')
    return b


def match(cell, build):
    """Key of the built-in annotation for a cell name and genome build, or None."""
    for key, e in STATES.items():
        if e['genome'] == build and re.match(rf'(?:{e["cell"]})(?:$|[_:.\-])', cell):
            return key
    return None


def resolve(spec, cell, build, refs):
    """(states path or None, states key or None, cytoband path or None, message) for one cell.
    spec: the `states` setting; refs: the folder that holds downloaded annotations."""
    if isinstance(spec, dict):
        spec = spec.get(cell)
    if not spec:
        return None, None, None, None
    band = CYTOBANDS.get(build)
    band_path = str(Path(refs) / 'annotations' / band['file']) if band else None
    if str(spec).lower() == 'auto':
        key = match(cell, build)
        if key is None:
            known = ', '.join(f"{k} ({e['genome']})" for k, e in STATES.items())
            return None, None, band_path, (f'states: auto: no built-in chromatin states for cell {cell!r} on genome '
                                           f'{build or "unknown"} (built in: {known}); give a BED file with '
                                           f'states: {{{cell}: FILE}}. Panels g-i are skipped.')
        return str(Path(refs) / 'annotations' / STATES[key]['file']), key, band_path, None
    return str(spec), None, None, None


def entry_of(filename):
    """The STATES or CYTOBANDS entry of a downloaded file name."""
    for e in list(STATES.values()) + list(CYTOBANDS.values()):
        if e['file'] == filename:
            return e
    raise KeyError(filename)


def fetch(filename, folder):
    """Download a built-in annotation file (its `file`, relative to folder) and check its md5; return its path.
    A file already there with the right md5 is kept. Uses download.py (resumes; md5 mismatch: an error)."""
    from . import download
    e = entry_of(filename)
    dest = Path(folder) / filename
    dest.parent.mkdir(parents=True, exist_ok=True)
    got = download.download(e['url'], e['md5'], connections=4, folder=dest.parent)
    assert Path(got) == dest, (got, dest)
    return dest


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_states(path, chrom):
    """[(start, end, state)] of one chromosome from a BED file, sorted by start."""
    rows = []
    with _open(path) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if len(x) >= 4 and x[0] == chrom:
                rows.append((int(x[1]), int(x[2]), x[3]))
    return sorted(rows)


def state_order(key, names):
    """(state, label) pairs in figure order: a built-in annotation's order, else the names in natural order
    (numeric prefixes as numbers: 1_..., 2_..., 10_...)."""
    if key:
        return [(s, lab) for s, lab in STATES[key]['order'] if s in names]
    nat = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]
    return [(s, s) for s in sorted(names, key=nat)]


def read_bands(path, chrom):
    """[(start, end)] of acen and gvar bands of one chromosome (UCSC cytoBand: chrom, start, end, name, stain)."""
    out = []
    with _open(path) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if len(x) >= 5 and x[0] == chrom and x[4] in ('acen', 'gvar'):
                out.append((int(x[1]), int(x[2])))
    return sorted(out)

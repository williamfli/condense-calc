"""Turn Park 2025's titration Source Data sheets into the pipeline's titration tables.

    python tools/titration_from_source_data.py OUTDIR XLSX [XLSX ...]

The sheets (Fig. 1b, ED Figs. 3d, 7c, 8a, 10a) are Park's titration files: blocks headed by
"<cell> <sample> <agent>", then a header row, then rows of
  conc, 3 replicate readings, 3 replicate fractions, mean fraction, std, Sample#
Rows whose Sample# is a number were sequenced; that number is the titration ID (tnum). Writes
OUTDIR/<cell>_<sample>_<agent>.tsv with columns tnum, conc, survival_fraction (the mean),
the names used in the GEO files (spermine -> sp, ...). Standard library only.
"""
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
AGENT = {'spermine': 'sp', 'spermidine': 'spd', 'cohex': 'CoH', 'peg': 'PEG', 'ca': 'Ca', 'mg': 'Mg',
         'hp1a': 'HP1a', 'hp1b/suv39h1': 'HP1bSUV'}
# Source Data names -> the names in GEO files and the pipeline's sample sheets
CELL = {'GM12878': 'GM', 'mCD8T (WT)': 'mCD8T:WT', 'mCD8T (DFMO)': 'mCD8T:DFMO', 'mCD8T (KO)': 'mCD8T:ODCKO'}
SAMPLE = {'reconstituted NCP': 'synNCP'}
# sheets whose block has no "<cell> <sample> <agent>" title row (from the figure legend)
UNTITLED = {'Extended Data Fig.3d': ('E14', 'NCP', 'spermine')}


def sheets(path):
    z = zipfile.ZipFile(path)
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        shared = [''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS)]
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    for s in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets', NS):
        target = rels[s.get(REL)].lstrip('/')
        target = target if target.startswith('xl/') else 'xl/' + target
        rows = []
        for row in ET.fromstring(z.read(target)).iter('{%s}row' % NS['m']):
            vals = []
            for c in row.findall('m:c', NS):
                v = c.find('m:v', NS)
                vals.append('' if v is None else (shared[int(v.text)] if c.get('t') == 's' else v.text))
            rows.append(vals)
        yield s.get('name'), rows


def blocks(rows):
    """Yield (cell, sample, agent, [(tnum, conc, frac)]) for each titration block."""
    i = 0
    while i < len(rows):
        r = [x.strip() for x in rows[i]]
        nxt = rows[i + 1] if i + 1 < len(rows) else []
        if len(r) >= 3 and r[2] and nxt and re.search(r'sample#', ' '.join(nxt), re.I):
            cell, sample, agent = r[0], r[1], r[2]
            out, i = [], i + 2
            while i < len(rows) and rows[i] and _num(rows[i][0]) is not None:
                v = rows[i]
                tnum = _num(v[-1])
                if tnum is not None and float(tnum).is_integer():
                    out.append((int(tnum), float(v[0]), float(v[7])))
                i += 1
            yield cell, sample, agent, out
        else:
            i += 1


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for path in sys.argv[2:]:
        for name, rows in sheets(path):
            if name in UNTITLED:
                rows = [list(UNTITLED[name])] + rows
            for cell, sample, agent, pts in blocks(rows):
                cell, sample = CELL.get(cell, cell), SAMPLE.get(sample, sample)
                if not pts:
                    continue
                ag = AGENT.get(agent.lower(), agent)
                f = out / f'{cell}_{sample}_{ag}.tsv'
                with open(f, 'w') as fh:
                    fh.write('tnum\tconc\tsurvival_fraction\n')
                    for t, c, s in sorted(pts):
                        fh.write(f'{t}\t{c:.6g}\t{s:.6g}\n')
                print(f'{Path(path).name} [{name}] -> {f.name}: tnums {[p[0] for p in sorted(pts)]}')


if __name__ == '__main__':
    main()

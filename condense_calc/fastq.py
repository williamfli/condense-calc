"""Cheap checks of FASTQ files before alignment: the first 10,000 reads of each file (milliseconds).

Settings `fastq_check`:
  lite  (default) stop if every base quality is '?': SRA Lite, which keeps no base qualities (ENA serves some
        runs so). bowtie2 aligns such reads worse, and input and supernatant libraries may differ in it.
  park  (Park mode) also stop unless the read names have the form Park et al. aligned, the instrument's:
        '@<spot name> <mate>:N:0:<index>' or '@<spot name> <mate>', the mate number one space after the
        name. bowtie2 picks among equally good hits by a seed from the read name, so ENA's
        '@SRR....<n> <spot name>/1' does not give Park's alignments.
  off   no check
Both point to the fix: FASTQs from NCBI's SRA Normalized data with instrument names (`condense fetch`).
"""
import gzip
import re

MODES = ('off', 'lite', 'park')
N_READS = 10_000
PARK_NAME = re.compile(r'@\S+ ([12])(:|$)')
FIX = ("Make the FASTQs from NCBI's SRA Normalized data: condense fetch SRR... (or fasterq-dump --split-files "
       "--seq-defline '@$sn $ri:N:0:$sg' --qual-defline '+')")
WHERE = 'in the --settings file of a sample-sheet run, or in the config file'


def _open(path):
    with open(path, 'rb') as f:
        gz = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gz else open(path)


def check_file(path, mate, mode='lite', n=N_READS):
    """Check the first n reads of one FASTQ file (mate 1 or 2). Returns what was seen, with a list of
    problems (empty if it passes)."""
    quals, first, bad_name, reads = set(), None, None, 0
    with _open(path) as f:
        for i, line in enumerate(f):
            if i >= 4 * n:
                break
            line = line.rstrip('\n')
            if i % 4 == 0:
                reads += 1
                first = first or line
                m = PARK_NAME.match(line)
                if bad_name is None and not (m and m.group(1) == str(mate)):
                    bad_name = line
            elif i % 4 == 3:
                quals.update(line)
    res = {'file': str(path), 'mate': mate, 'reads_checked': reads, 'first_name': first,
           'sra_lite': quals == {'?'}, 'park_names': bad_name is None, 'problems': []}
    if mode in ('lite', 'park') and res['sra_lite']:
        res['problems'].append(f"{path}: every base quality of the first {reads:,} reads is '?', i.e. SRA Lite, "
                               f'which has no base qualities. {FIX}. To align these reads anyway, set '
                               f'fastq_check: off ({WHERE}).')
    if mode == 'park' and bad_name is not None:
        res['problems'].append(f"{path}: read name {bad_name!r} is not in the form Park et al. aligned "
                               f"('@<spot name> {mate}:N:0:...', the mate number one space after the name); "
                               'bowtie2 places multi-mapping reads by their name, so the GEO tables are not '
                               f"reproduced. {FIX}. To align them anyway, set fastq_check: "
                               f"{'off' if res['sra_lite'] else 'lite'} ({WHERE}).")
    return res


def check(r1s, r2s, mode='lite', n=N_READS):
    """Check read-1 and read-2 files (lists of lanes). Returns one result per file (none with mode off)."""
    if mode not in MODES:
        raise ValueError(f'fastq_check must be one of {", ".join(MODES)}, not {mode!r}')
    if mode == 'off':
        return []
    return [check_file(f, mate, mode, n) for files, mate in ((r1s, 1), (r2s, 2)) for f in files]


def problems(results):
    return [p for r in results for p in r['problems']]

"""Paths for the step 0 write-up (scripts 03 and 04). Scripts 01 and 02 take theirs as arguments."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEP0 = HERE.parent                                   # geo_test/step0/
REPO = STEP0.parents[2]                               # ~/condense-seq
LAB = REPO / 'lab'
RESULTS, FIGURES = STEP0 / 'results', STEP0 / 'figures'
BIG = LAB / 'outputs/condense-seq-pipeline/geo_step0'  # large derived tables
GEO = LAB / 'data/GSE252941'
NETSCRATCH = Path('/n/netscratch/zhuang_lab/Lab/wfli/condense-seq')
RAW = NETSCRATCH / 'raw/GSE252941'
RUN = NETSCRATCH / 'pipeline_runs/geo_step0'          # the Snakemake run (bowtie2 seed 0)
SEED1 = RUN / 'seed1'                                 # 02: same reads, bowtie2 --seed 1
BT2_VERSIONS = ['2.3.5.1', '2.4.5']                   # 05a: same reads, older bowtie2 (RUN/bt2_<v>)
IDX4DN = RUN / 'idx_4dn'                              # 06: same reads, 4DN's GRCh38 bowtie2 index
MANIFEST = STEP0.parent / 'manifest.tsv'

COL = 'H1_NCP_PEG_6_2rep'                             # the one score column (ID 6 vs input)
KEY = ['Chromosome', 'Start', 'End']
CHROMS = [f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']]
# low-MAPQ classes: share of a bin's fragments (input + supernatant) with MAPQ < 10
LOWQ_EDGES = [-0.001, 0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
LOWQ_LABELS = ['0', '0-2%', '2-5%', '5-10%', '10-20%', '20-50%', '>50%']

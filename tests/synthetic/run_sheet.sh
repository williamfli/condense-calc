#!/bin/bash
# The same simulated experiment through a sample sheet (`condense SHEET`), compared with the
# config run of tests/synthetic/run.sh: every table must hold the same values (only names differ).
# Run from the repo root, in the condense-calc env, after run.sh OUTDIR [--park]:
#   bash tests/synthetic/run_sheet.sh OUTDIR [--park]       (the same mode as that run)
set -euo pipefail
OUT=$(realpath "$1"); shift
IN=$OUT/input
# sheet from the simulation's sample sheet and titration table: series SIM (shallow) and SIM_deep
awk -F'\t' 'NR==FNR { if (FNR > 1) { conc[$1] = $2; frac[$1] = $3 }; next }
            FNR == 1 { print "series\tpoint\tconc\tfraction\tfastq"; next }
            { s = ($6 == "deep") ? "SIM_deep" : "SIM"; print s "\t" $4 "\t" conc[$4] "\t" frac[$4] "\t" $7 }' \
  "$IN/titration/SIM_NCP_sp.tsv" "$IN/samples.tsv" > "$IN/sheet.tsv"
printf 'threads: 2\nbowtie2: bowtie2\n' > "$IN/settings.yaml"
python -m condense_calc "$IN/sheet.tsv" -o "$OUT/run_sheet" --genome "$IN/ref/sim.fa" --nucleosomes \
  --settings "$IN/settings.yaml" --cores "${SLURM_CPUS_PER_TASK:-2}" --quiet "$@"
python - "$OUT" <<'PY'
import gzip, sys
from pathlib import Path
out = Path(sys.argv[1])
names = {'SIM_NCP_sp_1rep_deep': 'SIM_deep', 'SIM_NCP_sp_1rep': 'SIM'}
body = lambda p: gzip.open(p, 'rt').read().split('\n', 1)[1]
bad = n = 0
for old in sorted((out / 'run' / 'tables').rglob('*.gtab.gz')):
    rel = str(old.relative_to(out / 'run' / 'tables'))
    key = next(k for k in names if k in rel)
    new = out / 'run_sheet' / 'tables' / rel.replace(key, names[key])
    same = new.exists() and body(new) == body(old)
    n += 1; bad += not same
    if not same:
        print('DIFFERENT', rel)
print(f'{n - bad}/{n} tables identical to the config run')
sys.exit(bad > 0 or n == 0)
PY

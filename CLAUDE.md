# condense-seq-pipeline

Repo-level context (`../CLAUDE.md`) and cluster rules (`~/CLAUDE.md`) load automatically. Read
`README.md` here first: stages, settings, validation, how to run on real data.

- **What it is:** raw condense-seq reads → GEO-format condensability tables. A Python 3 port of
  Park's `prepro_scripts`, byte-identical to the originals on synthetic data.
- **Faithful by default.** Any change to a stage's default behaviour must keep `tests/equivalence/`
  passing (identical to Park's scripts). New behaviour goes behind an opt-in option and into the
  README's "Departures" table.
- **After changing code**, run:
  - `pytest tests/test_units.py`
  - `sbatch tests/run_synthetic.sbatch`
  - `sbatch tests/equivalence/run_original.sbatch`; the comparison table must stay all identical.
- **Env:** `condseq` (`/n/holylabs/zhuang_lab/Lab/wfli/conda/envs/condseq`); spec in `env/`.
  Python 2.7 for the equivalence test comes from `module load python/2.7.16-fasrc01`.
- **Where things go:**
  - outputs: `../lab/outputs/condense-seq-pipeline/<run>/`
  - references: `../lab/data/ref/` (not installed yet)
  - logs: `../lab/logs/condseq_pipeline/`
- **Not here:** the Park polymer simulation (see `../aim1-crossmodal/`).

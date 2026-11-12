# CLAUDE.md

## What this project is
Pipeline for AlphaFold2 (ColabFold) models of putative nylon hydrolases (NylC homologs): prepare cleaved alpha/beta FASTA files,
post-process the models to match NylC X-ray structure 5y0m (chain IDs, residue numbering), and compute per-enzyme features for ML models.
Panel: NylC-GYAQ (reference, row 1 of the MSA) + nyl001-nyl095 (Zenodo 10.5281/zenodo.14976747) + dp001-dp189. The 24 extra picks nyl096-nyl119 are deliberately excluded.

## Layout
- `nylc_all_prep_20261111.ipynb`: MSA -> `results/cleavage_site.csv` + ColabFold FASTA files.
- `nylc_all_postprocess_20261111.ipynb`: renumber models -> features -> `results/nylc_all.merged_data_<date>.csv` -> comparison with `legacy/`.
- `scripts/`: PyMOL post-processing, BSA, fpocket (called by the notebook; PyMOL runs in subprocesses because it can kill the kernel).
- `data/`: small inputs (MSA, sequences, activity, 5y0m). Never write here. `legacy/`: reference outputs for regression checks.
- `af2_models/raw/`: relaxed rank-1 PDB + `scores` JSON files from ColabFold (large, gitignored, download separately; see README).
- `results/`: all outputs (gitignored except what you choose to commit).

## Environment
`mamba env create -f environment.yml && mamba activate nylc`. Run notebooks from the repo root. Paths are set via
`NYLC_DATA_DIR`, `NYLC_RESULTS_DIR`, `NYLC_RAW_AF2_DIR`, `NYLC_LEGACY_DIR`, `NYLC_N_JOBS`, `NYLC_FORCE=1` (recompute).
Headless: `jupyter nbconvert --to notebook --execute <nb> --output-dir results/executed --ExecutePreprocessor.timeout=-1`.

## Conventions and gotchas
- Enzyme names are lowercase, three-digit (`nyl001`, `dp001`); old files may use `nyl01`. Match names exactly, never by substring/glob prefix.
- Residue numbers: alpha chain 1..n, beta chain continues at `Thr267` column of `cleavage_site.csv` (= alpha length + 1). Chains A-D after processing (alpha+beta share an ID).
- MSA features (pLDDT regions, loops, key residues) use alignment columns. Do not re-align or drop all-gap columns when removing sequences, or features change.
- Raw ColabFold chains: A-D alpha, E-H beta. The post-processing script asserts this.
- Never modify the original notebooks or the legacy files in `nylc_all/data/`. Skip-flag files (`*_DONE`) were replaced by output-exists checks.
- Legacy column `Thr219_equiv` is `Thr267_equiv` here. See `AUDIT_20261111.md` for all behavioral differences from the legacy notebooks.

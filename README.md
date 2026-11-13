# NylC homolog AlphaFold2 models: post-processing and feature calculation

Pipeline to prepare sequences of putative nylon hydrolases (NylC homologs) for AlphaFold2/ColabFold, post-process the models
(canonical chain IDs and residue numbering matching the NylC X-ray structure, PDB 5y0m) and compute features for machine learning.

Panel: NylC-GYAQ (reference) + 95 homologs of the original Diversity Panel from Drufva et al. Chem. Catal. 2025. DOI: 10.1016/j.checat.2025.101418 (`nyl001`-`nyl095`) plus 189 additional homologs (`dp001`-`dp189`) = 285 models.

Authors: Megan M. Wisniewski and Jerry M. Parks

## Install
```bash
mamba env create -f environment.yml
mamba activate nylc
```
(Linux or macOS; PyMOL, fpocket and MAFFT come from conda-forge.)

## Get the AlphaFold2 models (not stored in this repo)
Raw ColabFold output (relaxed rank-1 PDB + `scores` JSON files) for all models is on Zenodo, DOI [10.5281/zenodo.24201789](https://zenodo.org/records/24201789):

| File | Contents | MD5 |
|---|---|---|
| `raw_af2_dp95_zenodo.tgz` | NylC-GYAQ and 95 homologs (`nyl001`-`nyl095`) | `e07213ae8a29967e6deaab36308fb30a` |
| `raw_af2_dp189_zenodo.tgz` | 189 additional homologs (`dp001`-`dp189`) | `19528067d81e532d9e286a87689b89d9` |

```bash
mkdir -p af2_models/raw
for f in raw_af2_dp95_zenodo raw_af2_dp189_zenodo; do
  curl -L -o $f.tgz "https://zenodo.org/records/24201789/files/$f.tgz?download=1"
  tar xzf $f.tgz && mv $f/* af2_models/raw/ && rm -r $f $f.tgz
done
```
The notebooks expect files named `<name>_4x2_relaxed_rank_001_*.pdb` and `<name>_4x2_scores_rank_001_*.json` in `af2_models/raw/` (override with `NYLC_RAW_AF2_DIR`).
The notebooks generate the processed models (renumbered, chain-renamed) for all 285 enzymes from these raw files.

## Run
1. `nylc_all_prep_20261111.ipynb` writes `results/cleavage_site.csv` and the ColabFold FASTA files (`results/fasta_af2/`). Run ColabFold on those if you want to regenerate models.
2. `nylc_all_postprocess_20261111.ipynb` writes processed models, per-feature CSVs (`results/features/`) and `results/nylc_all.merged_data_<date>.csv`; it ends with a comparison against `legacy/nylc_all.merged_data_20250310.csv`.

Headless:
```bash
jupyter nbconvert --to notebook --execute nylc_all_prep_20261111.ipynb --output-dir results/executed
jupyter nbconvert --to notebook --execute nylc_all_postprocess_20261111.ipynb --output-dir results/executed --ExecutePreprocessor.timeout=-1
```
Runtime is about 10 minutes with 24 workers (`NYLC_N_JOBS`).

## Repository contents
| Path | Purpose |
|---|---|
| `data/` | MSA (`nylc_all.mafft.fasta`, L-INS-i via Jalview), sequences, activity table, 5y0m tetramer |
| `scripts/` | PyMOL post-processing, BSA, fpocket descriptors, derivation of `data/` from the legacy files |
| `legacy/` | Legacy merged feature table and cleavage sites used for regression checks |
| `AUDIT_20261111.md` | Errors found in the original notebooks and how they were handled |

## Notes
- BSA and fpocket descriptors depend slightly on the coordinate frame. The processing script therefore rotates all models into a fixed frame; do not change it if you want to reproduce published values.
- Legacy column `Thr219_equiv` is named `Thr267_equiv` here.

# cancer_drug_response_prediction_multiomics

Predicting cancer drug response from multi-omics data (gene expression, copy number, mutations):
train on cancer cell lines (CCLE, GDSC), predict patient-derived xenografts (PDX).
Based on the compgen2021 capstone project (MDC/BIMSB, Akalin lab), redone from scratch.

## Links
- Task description: https://bimsbstatic.mdc-berlin.de/akalin/CompGenom2021/compgen2021Capstone/doc/CompGen2021CapstoneProject.html
  (local copy: `../compgen2021/capstone/instructions/task_description/CompGen2021CapstoneProject.html`)
- Capstone repo: https://github.com/BIMSBbioinfo/compgen2021Capstone
- Data archive: https://bimsbstatic.mdc-berlin.de/akalin/AAkalin_DrugResponse/prepared/AAkalin_DrugResponse-prepared.tar.gz

## Tasks
1. **Multi-omics:** expression + CNV + mutations (`all_intersect_genes`). Train on each cell-line dataset, test on the others, then predict PDX.
2. **Gene panel:** CNV + mutations of the 882 oncoKB genes only (`oncoKB_intersect_genes`).
3. **Batch correction:** remove technical differences between cell lines and PDX, then compare PDX predictions with Tasks 1 and 2.

## Layout
```
cancer_drug_response_prediction_multiomics/
├── data/
│   ├── raw/          copy of the original capstone data, never modified (not in git)
│   └── processed/    cleaned matrices written by scripts/ (not in git)
├── R/                reusable functions, sourced by the scripts
├── scripts/          analysis steps, numbered in run order
├── results/
│   ├── figures/
│   └── tables/
└── models/           fitted models (not in git)
```
`data/raw/` is a copy of `../compgen2021/capstone/data/`. The task docs, the 2021 report and the
rebuilt 2021 scripts stay in `../compgen2021/capstone/`.
Run everything from the project root, so paths look like `data/raw/CCLE/...`.
R environment: mamba `Env1` (R 4.5.1).

## Data (`data/raw/`)
| Dataset | Samples (all 3 omics) | Role |
|---|---|---|
| CCLE | 1,291 | training |
| GDSC_broad | 507 | training / testing |
| GDSC_sanger | 281 | training / testing |
| PDX_all | 399 | **testing only**: never train or tune on it |

Each has `all_intersect_genes/` (gex, cnv, mut), `oncoKB_intersect_genes/` (cnv, mut) and
`drug_response/` (one file per drug: ALPELISIB, BUPARLISIB, CETUXIMAB, ERLOTINIB, RIBOCICLIB, TRAMETINIB).
Response values are AUC: lower = more sensitive. Format details: `data/raw/README`.

## Known data pitfalls
1. **PDX expression is untransformed** (max ≈ 325,234); cell-line expression is already log-transformed (0 to ≈ 16).
2. **PDX drug files are mostly empty**: Erlotinib has 29 real values, Trametinib 39.
3. **Ribociclib** has only 42 / 36 / 3 cell lines (CCLE / GDSC_broad / GDSC_sanger), with nearly constant values.
4. **GDSC_broad and GDSC_sanger drug files are identical**; the datasets differ only in omics and cell lines.
5. **CCLE and GDSC_broad share cell lines** (544 of 545), so testing one on the other is not independent.
6. **GDSC has duplicate measurements** for some cell lines (Erlotinib, Trametinib).
7. **Mutation formats differ**: CCLE has `Variant_Classification`; it is empty in GDSC; PDX writes `C1277G`
   instead of `p.C1277G` and has no gene names. Gene-level 0/1 encoding avoids most mismatches.

## Not included
`ccle_mut.rds` used in the course's Simple Example PDF was shared separately and is not in the public archive.

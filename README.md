# Cancer drug response prediction from multi-omics data

Can drug response measured in cancer cell lines predict how real tumours respond?
This project trains machine-learning models on three cell line drug screens and tests them on
patient-derived xenografts (PDX: human tumours grown in mice), for six targeted cancer drugs:

- **CCLE**: Cancer Cell Line Encyclopedia (Broad Institute);
- **GDSC-Broad** and **GDSC-Sanger**: Genomics of Drug Sensitivity in Cancer, cell lines profiled
  at the Broad and at the Sanger Institute;
- **PDX**: Novartis PDX Encyclopedia, used only for testing.

The datasets differ systematically (batch effect), so they are batch-corrected first. Two feature
sets are compared: full multi-omics (expression, copy number, mutations of ~18,000 genes) and a
clinical-style **gene panel** (copy number and mutations of 882 cancer genes, no expression).

PDX is used only as a test set: no model, setting or ensemble is trained or chosen on it. The
whole analysis is a Snakemake pipeline (see How to run).

## Main results

| | Multi-omics | Gene panel |
|---|---|---|
| Cross-validation R², mean of best models | **0.21** | 0.10 |
| R² on another cell line dataset, mean | **0.27** | 0.14 |
| PDX: Spearman with response (ensemble, mean over drugs) | **0.17** | 0.05 |
| PDX: ROC AUC, tumour kept growing vs not (ensemble, mean) | **0.62** | 0.53 |

1. **Batch correction removes most of the dataset differences** between the cell line datasets
   (Figure 1a–c). The results hardly depend on it:
   - Without correction, cross-validation is identical.
   - Transfer between cell line datasets is slightly lower (R² 0.25 vs 0.27).
   - PDX results are unchanged.

   `results/tables/summary.csv` has both.
2. **Response is partly predictable within a dataset.** Cross-validated R² ranges from ~0.1
   (Buparlisib, Cetuximab) to ~0.5 (Trametinib). SVM, elastic net and PLS did best of six model
   types (Figure 1d–g).
3. **Expression matters.** The gene panel halves performance everywhere and is near chance on PDX
   (Figures 1f–g, 2, 3).
4. **On PDX, the models carry signal for some drugs only.** Tumours that kept growing under
   Erlotinib or Trametinib get clearly higher predicted resistance (ROC AUC 0.86 and 0.79).
   Buparlisib, Cetuximab and Ribociclib are at chance. Erlotinib (25 PDX) and Trametinib (37)
   have few samples, so these numbers are uncertain (Figure 3).
5. **The models pick up known biology** (Figure 4):
   - *PIK3CA* mutation predicts sensitivity to the PI3Kα inhibitor Alpelisib.
   - *EGFR* expression and its ligand *BTC* predict sensitivity to Erlotinib.
   - *SPRY2*, a readout of MAPK pathway activity, predicts sensitivity to the MEK inhibitor
     Trametinib.
6. **R² on PDX is always strongly negative (about −9)** because the PDX response is a tumour
   volume change, not a dose-response AUC, so the values sit on a different scale. Correlations
   and ROC AUC are the fair measures there.

## Figures

All figures are drawn by `scripts/08_figures.R` from the tables in `results/tables/`. All
results are from batch-corrected data. The exception is the model comparison (Figure 1d, e),
which gives the same result either way. Figure 1a and c also show the uncorrected data.

![Figure 1](results/figures/fig1_data_and_cv.png)

**Figure 1. Batch effect between datasets, model choice and cross-validated performance.**
- (**a**, **b**) First two principal components of gene expression for all 2,478 samples:
  1,291 CCLE, 507 GDSC-Broad and 281 GDSC-Sanger cell lines and 399 PDX tumours. The 2,000 genes
  with the highest variance over all samples are used, scaled, the same in both panels.
  Percentages are the variance explained.
  - (**a**) Uncorrected, the datasets separate, although all GDSC-Broad cell lines are also in
    CCLE.
  - (**b**) After batch correction (`batchelor::regressBatches`, which removes each dataset's own
    mean per gene) the cell line datasets overlap. PDX stays partly apart, which is expected
    biologically: they are tumours grown in mice, while the others are cells grown in plastic.
- (**c**) Mixing: the share of each sample's 20 nearest neighbours (in the first 10 PCs) that
  come from its own dataset. Grey: uncorrected; black: batch-corrected. Dashed line: the 35%
  expected if the datasets were fully mixed. For cell lines it falls from 85% to 55%,
  for PDX from 99% to 69%.
- (**d**) Six model types compared by 10-fold cross-validated *R*² (1 − SSE/SST on out-of-fold
  predictions), one point per drug × cell line dataset (17 combinations; Ribociclib in
  GDSC-Sanger has too few samples). Multi-omics features, the same folds for every model.
  Black bars: mean.
  - This comparison was run before batch correction. Within one dataset the correction only
    shifts each feature by a constant, which does not change cross-validation.
  - Kept for all later steps (blue): SVM with radial kernel (0.19), elastic net (0.18) and PLS
    regression (0.17).
  - Dropped (grey): random forest (0.17), gradient boosting (xgboost, 0.15) and k-nearest
    neighbours (0.14).
- (**e**) Number of the 17 combinations in which each model type had the highest *R*².
- (**f**, **g**) Cross-validated *R*² of the best of the three kept model types per drug and
  dataset, with (**f**) multi-omics features and (**g**) the gene panel.
  - The model type is shown under each value: EN, elastic net; SVM, support vector machine;
    PLS, partial least squares. n.t.: not trained (3 samples).
  - The panel roughly halves performance (mean *R*² 0.21 vs 0.10).

![Figure 2](results/figures/fig2_cross_dataset.png)

**Figure 2. Transfer of models between cell line datasets.**
- (**a**) *R*² when the best model of one dataset predicts the other two.
  - Points: one train → test pair (up to six per drug, two for Ribociclib; test sets with ≥ 10
    samples only).
  - Bars: mean.
- (**b**, **c**) Mean *R*² over drugs for each training → test pair, with (**b**) multi-omics
  features and (**c**) the gene panel.
  - CCLE and GDSC-Broad share their cell lines and nearly the same drug responses, so transfer
    between them is optimistic.
  - Tests involving GDSC-Sanger are the more independent ones (54 cell lines shared with CCLE).
  - Models trained on GDSC-Sanger, the smallest dataset, transfer worst.

![Figure 3](results/figures/fig3_pdx.png)

**Figure 3. Prediction of PDX tumour response with models trained on cell lines.**
For each drug, the PDX response is predicted with the ensemble of the best models of the three
cell line datasets. No model, setting or ensemble member was chosen using PDX.
- (**a**) Spearman correlation between predicted cell line AUC and the PDX response; n is the
  number of PDX models with a response.
- (**b**) ROC AUC for separating tumours that kept growing under treatment (mRECIST progressive
  disease, PD) from the rest; 0.5 is chance. Response classes come from Xeva (see Data).
  PD / all classified models: Alpelisib 48/137, Buparlisib 53/164, Cetuximab 38/62,
  Erlotinib 8/20, Ribociclib 71/164, Trametinib 9/36.
- (**c**) Predicted AUC of the multi-omics ensemble by response class.
  Boxes: median and interquartile range; points: PDX models.

The PDX response is a change in tumour volume, not a dose-response AUC, so *R*² is not
meaningful here (about −9 everywhere; see `pdx_metrics.csv`). Only Erlotinib and Trametinib
are clearly above chance, and both have few PDX models.

![Figure 4](results/figures/fig4_top_features.png)

**Figure 4. Strongest features of the elastic net models.**
Ten features with the largest mean standardized weight per drug, from the elastic net models
(multi-omics) of each cell line dataset. Elastic net is the one kept model type with
a weight per feature that can be read directly.
- **Standardized weight:** the model coefficient × the feature's standard deviation in the
  training samples, i.e. the change in predicted AUC for a one-SD change in the feature. It is
  averaged over the datasets, with 0 where a model did not select the feature.
- **Colours:** blue, a higher value predicts sensitivity (lower AUC); red, it predicts
  resistance.
- **Numbers:** how many of the drug's models selected the feature. Features selected in every
  dataset are the more reliable markers.
- **Labels:** expr, expression; CNV, copy number; mut, mutation. Gene symbols come from
  org.Hs.eg.db.

Full list (top 20 per drug): `results/tables/top_features.csv`.

## Methods

The steps are rules in the `Snakefile`, one script each:

| Step | Script | What it does |
|---|---|---|
| 1 | `01_prepare_data.R` | Reads all datasets; log2(x + 1)-transforms PDX expression (raw FPKM, unlike the cell lines); encodes mutations as gene-level 0/1 (non-silent); averages repeated drug measurements; saves multi-omics and panel matrices and the Xeva PDX response measures |
| 2 | `02_batch_correction.R` | `batchelor::regressBatches` per omic and feature set: removes each dataset's own mean per feature. Uses PDX omics, never PDX responses. Also the PCA before / after |
| 3 | `03_model_comparison.R` | Six model types (elastic net, random forest, kNN, SVM, PLS, xgboost), 10-fold CV with the same folds, multi-omics data. Keeps the top 3 |
| 4 | `04_train_models.R` | One model per drug × cell line dataset × feature set, 10-fold CV, best of the 3 types kept. Also trained on uncorrected data, for the comparison in `summary.csv` |
| 5 | `05_cross_dataset.R` | Each best model applied to the other two cell line datasets |
| 6 | `06_predict_pdx.R` | PDX prediction with a single best model and an ensemble |
| 7 | `07_top_features.R` | Standardized weights of the elastic net models |
| 8 | `08_figures.R` | All figures (from the result tables) and `results/tables/summary.csv` (all results, with and without batch correction) |

Details:
- **Features**, chosen on the training samples only: top 5000 most variable expression genes, top
  5000 most variable CNV genes, genes mutated in ≥ 5% of samples. Gene panel: all CNV genes and
  mutations down to 1% (curated cancer genes, e.g. EGFR mutations for Erlotinib).
- **R²** is 1 − SSE/SST on out-of-fold predictions (not caret's squared correlation), so it can be negative.
- **Single best model** per drug: of the best model of each cell line dataset, the one with the
  highest R² on the other cell line datasets. **Ensemble**: average of those models, leaving out
  any with CV R² ≤ 0.
- **PDX evaluation**: R², Pearson and Spearman against the PDX response; and ROC AUC for
  separating tumours that kept growing (mRECIST "progressive disease") from the rest, using
  [Xeva](https://bioconductor.org/packages/Xeva/) response classes.
- **Batch correction**: per-dataset z-scoring and fastMNN were tried too and did no better
  (fastMNN mixed the datasets most but predicted PDX worst). regressBatches corrects each dataset
  with its own samples only, so new data can be corrected the same way without retraining.
- Models with fewer than 20 training samples are skipped, test sets with fewer than 10 too.

## Data

Original data: compgen2021 capstone project (MDC/BIMSB, Akalin lab),
[archive](https://bimsbstatic.mdc-berlin.de/akalin/AAkalin_DrugResponse/prepared/AAkalin_DrugResponse-prepared.tar.gz),
[task description](https://bimsbstatic.mdc-berlin.de/akalin/CompGenom2021/compgen2021Capstone/doc/CompGen2021CapstoneProject.html).

| Dataset | Folder | What it is | Samples (multi-omics) | Role |
|---|---|---|---|---|
| CCLE | `CCLE/` | Cancer Cell Line Encyclopedia: a large panel of cancer cell lines, molecularly profiled and drug-screened at the Broad Institute | 1,291 (1,722 with panel data) | training / testing |
| GDSC-Broad | `GDSC_broad/` | Genomics of Drug Sensitivity in Cancer (Sanger Institute) drug screens, for the cell lines whose omics come from the Broad Institute. The cell lines are all in CCLE too, with separately processed omics and nearly the same drug responses (r = 0.96–1.00 for five drugs) | 507 | training / testing |
| GDSC-Sanger | `GDSC_sanger/` | GDSC drug screens for the cell lines profiled at the Sanger Institute; no overlap with GDSC-Broad, 54 cell lines also in CCLE | 281 | training / testing |
| PDX | `PDX_all/` | Patient-derived xenografts from the Novartis PDX Encyclopedia (Gao et al., *Nat Med* 2015): pieces of patient tumours grown and treated in mice. Closer to real tumours than cell lines grown in plastic | 399 | **testing only** |

Drugs: Alpelisib, Buparlisib (PI3K), Cetuximab, Erlotinib (EGFR), Ribociclib (CDK4/6), Trametinib (MEK).
Cell line response = AUC of the dose-response curve; lower = more sensitive. The PDX response
is a change in tumour volume (see Xeva below).

| Folder in `data/raw/` | What it is | Used |
|---|---|---|
| `CCLE/`, `GDSC_broad/`, `GDSC_sanger/`, `PDX_all/` | Per dataset: `all_intersect_genes/` (expression, CNV, mutations of the ~18,000 genes shared by all datasets), `oncoKB_intersect_genes/` (CNV and mutations of the 882 OncoKB cancer genes, the gene panel) and `drug_response/` | yes |
| `auxiliary/` | Gene lists from the course: `all_intersect_genes.tsv` (18,355 Ensembl IDs shared by all datasets), `oncoKB_intersect_genes.tsv` (the 882 panel genes) and `oncoKB_gene_list.tsv` (all 1,064 OncoKB genes with symbol, Entrez and Ensembl IDs) | no: the omics files are already limited to these genes |
| `xeva/` | `Xeva_PDXE_tsv.tar.gz`: the Novartis PDX Encyclopedia curated with the Bioconductor package [Xeva](https://bioconductor.org/packages/Xeva/) (ORCESTRA, [Zenodo record 19354439](https://zenodo.org/records/19354439)) | yes, for the PDX response classes |

**Xeva.** The course gives one PDX response value per tumour and drug. Xeva has the same PDX
models with the tumour volume curves and several response measures. `01_prepare_data.R` reads two
of its tables: `model.tsv` (model → patient and drug) and `sensitivity_model.tsv` (response per
model). It saves the mRECIST class, best response, best average response, slope and AUC of the
growth curve for the 399 PDX samples and six drugs (`data/processed/PDX_xeva_response.rds`).
They serve two purposes:
- They show what the course PDX value is: Xeva's best % change in tumour volume, rescaled
  (Spearman 0.93–0.99). It cannot be on the same scale as the cell line AUC.
- They give the mRECIST classes used in Figure 3b, c: CR / PR (tumour shrank), SD (stable) and
  PD (kept growing).

Xeva also has its own PDX omics, including continuous copy numbers. The course PDX copy numbers
are almost all neutral (~1% gains or losses, vs 40–55% in cell lines). Rebuilding them from Xeva
did not improve PDX prediction, so the models use only the course omics.

## How to run

Needs R 4.5 and Snakemake (≥ 8), with these R packages:
- **CRAN:** data.table, caret, glmnet, kernlab, pls, ranger, xgboost, doParallel, ggplot2,
  patchwork, ragg.
- **Bioconductor:** batchelor, BiocNeighbors, SummarizedExperiment, AnnotationDbi, org.Hs.eg.db.

1. Put the course data in `data/raw/` (`CCLE/`, `GDSC_broad/`, `GDSC_sanger/`, `PDX_all/`) and
   `Xeva_PDXE_tsv.tar.gz` in `data/raw/xeva/`.
2. From the project root, with R and Snakemake on the PATH:
   ```
   snakemake -n            # show what would run
   snakemake --cores 12    # run everything that is missing or out of date
   ```
   The model comparison and the training take about an hour each on 12 cores; the rest take
   minutes. A step reruns only when its script, its inputs or the R/ helpers it uses change, so
   editing `scripts/08_figures.R` only redraws the figures. Logs: `results/logs/`.

```
├── Snakefile          the pipeline
├── data/raw/          original data, never modified (not in git): 4 datasets, auxiliary/, xeva/
├── data/processed/    cleaned and batch-corrected matrices (not in git)
├── R/                 functions used by the scripts
├── scripts/           analysis steps, in run order
├── models/            fitted models, one folder per feature set and correction (not in git)
└── results/           figures/, tables/, logs/
```

## Limitations
- Cell line AUC and PDX tumour response measure different things, so PDX results are judged on
  ranking only.
- PDX sample sizes are small for Erlotinib and Trametinib, and responders (shrinking tumours) are
  rare for every drug (3–12 per drug), so PD vs not PD is used instead.
- The gene panel has more CCLE cell lines than the multi-omics data (no expression needed),
  which slightly favours it.
- The feature filter is fitted once per training set, not inside each CV fold, and the best model
  type is chosen by CV R²; both make CV R² slightly optimistic. Cross-dataset and PDX results are not affected.
- The top features come from one linear model type and show association, not causation.

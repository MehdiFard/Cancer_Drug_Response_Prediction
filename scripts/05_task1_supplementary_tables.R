# Task 1: collect all result tables into one Excel file, like a supplementary file.
# First sheet explains what each table is, second sheet explains every column.
#
# Run from the project root, after 02, 03 and 04.

library(data.table)
library(openxlsx)

out_file <- "results/tables/task1_supplementary_tables.xlsx"

# ---- What each table is ----
# sheet name (max 31 characters), csv file, description
tables <- rbindlist(list(
  list("S1_cv_method_summary", "task1_cv_method_summary.csv",
       "Task 1.1. Six model types compared with 10-fold CV on the same folds, averaged over all drug x dataset combinations. The top 3 by mean CV R-squared were kept."),
  list("S2_cv_metrics", "task1_cv_metrics.csv",
       "Task 1.1. 10-fold CV performance of every model type for every drug and cell line dataset, computed on the out-of-fold predictions."),
  list("S3_cv_r2_best", "task1_cv_r2_best.csv",
       "Task 1.1. CV R-squared of the best model (among the kept model types) per drug and dataset. Same values as the heatmap figure. Empty = not trained (fewer than 20 samples)."),
  list("S4_cv_predictions", "task1_cv_predictions.csv",
       "Task 1.1. Out-of-fold predictions: each sample predicted by the model of the CV fold that did not see it. One row per sample, drug, dataset and model type."),
  list("S5_cross_dataset", "task1_cross_dataset_metrics.csv",
       "Task 1.2. The best model of each cell line dataset applied to the other two datasets, for each drug."),
  list("S6_cross_dataset_summary", "task1_cross_dataset_summary.csv",
       "Task 1.2. Mean and median of the cross-dataset R-squared values per drug (up to 6 values: 3 training x 2 test datasets)."),
  list("S7_cross_dataset_by_dataset", "task1_cross_dataset_by_dataset.csv",
       "Task 1.2. Cross-dataset R-squared averaged by training dataset and by test dataset, to see if one dataset gives poor models or is hard to predict."),
  list("S8_pdx_model_choice", "task1_pdx_model_choice.csv",
       "Task 1.3. Candidate models for PDX prediction (the best model of each cell line dataset) and how the single best model and the ensemble members were chosen. Chosen from cell line results only, never from PDX."),
  list("S9_pdx_metrics", "task1_pdx_metrics.csv",
       "Task 1.3. Performance on PDX of the single best model and of the ensemble (average of the candidate models) per drug."),
  list("S10_pdx_predictions", "task1_pdx_predictions.csv",
       "Task 1.3. PDX predictions per sample and drug: observed response, single best model, ensemble, and each candidate model.")
))
setnames(tables, c("sheet", "file", "description"))

# ---- What each column is ----
# Columns that appear in several tables are described once in `common`
common <- c(
  drug        = "Drug name.",
  dataset     = "Cell line dataset (CCLE, GDSC_broad or GDSC_sanger).",
  method      = "Model type: glmnet = elastic net, ranger = random forest, knn = k-nearest neighbours, svmRadial = support vector machine with radial kernel, pls = partial least squares, xgboost = gradient boosted trees.",
  n           = "Number of samples the metrics are computed on.",
  r2          = "R-squared = 1 - SSE/SST (coefficient of determination). 1 = perfect, 0 = as good as predicting the mean, below 0 = worse than the mean. Not caret's default (squared correlation).",
  pearson     = "Pearson correlation between observed and predicted response.",
  spearman    = "Spearman (rank) correlation between observed and predicted response.",
  rmse        = "Root mean squared error, in AUC units.",
  sample_id   = "Sample (cell line or PDX model) ID.",
  obs         = "Observed drug response (AUC; lower = more sensitive).",
  pred        = "Predicted drug response (AUC).",
  mean_r2     = "Mean R-squared.",
  median_r2   = "Median R-squared.",
  CCLE        = "Value for the CCLE dataset / from the CCLE model.",
  GDSC_broad  = "Value for the GDSC_broad dataset / from the GDSC_broad model.",
  GDSC_sanger = "Value for the GDSC_sanger dataset / from the GDSC_sanger model."
)

# Columns specific to one table (or with a different meaning there)
specific <- list(
  S1_cv_method_summary = c(
    mean_rank = "Mean rank of the model type within each drug x dataset (1 = best of the six).",
    n_best    = "Number of drug x dataset combinations where this model type ranked first.",
    kept      = "TRUE if the model type is in the top 3 by mean CV R-squared and used in the next steps."),
  S2_cv_metrics = c(
    tune = "Tuning parameters chosen by CV (lowest RMSE).",
    rank = "Rank of this model type within the drug x dataset (1 = highest CV R-squared).",
    kept = "TRUE if the model type is one of the 3 kept ones.",
    best = "TRUE for the best kept model of this drug and dataset, used in Tasks 1.2 and 1.3."),
  S3_cv_r2_best = c(
    CCLE        = "CV R-squared of the best model trained on CCLE.",
    GDSC_broad  = "CV R-squared of the best model trained on GDSC_broad.",
    GDSC_sanger = "CV R-squared of the best model trained on GDSC_sanger."),
  S4_cv_predictions = c(
    pred = "Out-of-fold predicted response (AUC), from the best tuning parameters."),
  S5_cross_dataset = c(
    train    = "Dataset the model was trained on.",
    test     = "Dataset the model was applied to.",
    method   = "Model type of the best model of the training dataset.",
    n_shared = "Number of test samples that were also in the training set (same cell line in both datasets). High for CCLE vs GDSC_broad (Pitfall 5), so those tests are not independent."),
  S6_cross_dataset_summary = c(
    n_pairs      = "Number of train -> test pairs (6, or 2 for Ribociclib since GDSC_sanger has too few samples).",
    mean_pearson = "Mean Pearson correlation over the pairs."),
  S7_cross_dataset_by_dataset = c(
    dataset = "Cell line dataset.",
    role    = "train = averaged over the models trained on this dataset; test = averaged over the predictions made on this dataset."),
  S8_pdx_model_choice = c(
    dataset     = "Dataset the candidate model was trained on.",
    method      = "Model type of the candidate model.",
    cv_r2       = "CV R-squared of the candidate model (Task 1.1).",
    cross_r2    = "Mean R-squared of the candidate model on the other cell line datasets (Task 1.2).",
    single      = "TRUE for the single best model: the candidate with the highest cross_r2.",
    in_ensemble = "TRUE if the candidate is in the ensemble (CV R-squared > 0)."),
  S9_pdx_metrics = c(
    model  = "single = single best model; ensemble = average of the ensemble members.",
    models = "Which models (training dataset and model type) were used.",
    r2     = "R-squared = 1 - SSE/SST on PDX. Strongly negative because PDX AUC (from tumor volume) is on another scale than cell line AUC; use the correlations instead."),
  S10_pdx_predictions = c(
    obs         = "Observed PDX response (AUC from tumor volume; lower = more sensitive).",
    single      = "Prediction of the single best model.",
    ensemble    = "Prediction of the ensemble (average of the ensemble members).",
    CCLE        = "Prediction of the best CCLE model.",
    GDSC_broad  = "Prediction of the best GDSC_broad model.",
    GDSC_sanger = "Prediction of the best GDSC_sanger model (empty for Ribociclib, no model).")
)

# ---- Build the workbook ----
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
head_style  <- createStyle(textDecoration = "bold", fgFill = "#dde8f6", border = "bottom", wrapText = TRUE, valign = "top")
wrap_style  <- createStyle(wrapText = TRUE, valign = "top")
num_style   <- createStyle(numFmt = "0.000")
title_style <- createStyle(fontSize = 14, textDecoration = "bold")

dictionary <- list()
data_list  <- list()
for (i in seq_len(nrow(tables))) {
  sh  <- tables$sheet[i]
  dat <- fread(file.path("results/tables", tables$file[i]))
  data_list[[sh]] <- dat

  desc <- sapply(names(dat), function(col) {
    if (!is.null(specific[[sh]]) && col %in% names(specific[[sh]])) specific[[sh]][[col]] else common[col]
  })
  if (anyNA(desc)) stop("No description for column(s) ", paste(names(dat)[is.na(desc)], collapse = ", "), " in ", sh)
  dictionary[[sh]] <- data.table(sheet = sh, column = names(dat), description = unname(desc))
}
dictionary <- rbindlist(dictionary)

# README sheet
addWorksheet(wb, "README")
readme_top <- c(
  "Supplementary tables, Task 1: predicting drug response from multi-omics data",
  "Project: cancer_drug_response_prediction_multiomics (compgen2021 capstone, Task 1).",
  "Data: gene expression, copy number and mutations (all_intersect_genes) of three cell line datasets (CCLE, GDSC_broad, GDSC_sanger) and PDX; response = AUC per drug, lower = more sensitive.",
  "Models were trained on cell lines only (10-fold CV, one model per drug and dataset). PDX was used only as a test set, never to train or choose anything.",
  "Features per model: top 5000 most variable expression genes, top 5000 most variable CNV genes, genes mutated in >= 5% of samples; selected on the training samples only.",
  "R-squared is 1 - SSE/SST throughout (not caret's squared correlation).",
  "Sheet 'Columns' describes every column of every table.",
  paste("Generated:", format(Sys.Date()))
)
writeData(wb, "README", readme_top, startRow = 1, colNames = FALSE)
addStyle(wb, "README", title_style, rows = 1, cols = 1)
index <- rbind(tables[, .(sheet, description, source_file = file.path("results/tables", file),
                          rows = sapply(data_list[sheet], nrow))])
writeData(wb, "README", index, startRow = length(readme_top) + 2, headerStyle = head_style)
addStyle(wb, "README", wrap_style, rows = length(readme_top) + 2 + seq_len(nrow(index)), cols = 2, gridExpand = TRUE)
setColWidths(wb, "README", cols = 1:4, widths = c(30, 90, 45, 8))

# Columns sheet
addWorksheet(wb, "Columns")
writeData(wb, "Columns", dictionary, headerStyle = head_style)
addStyle(wb, "Columns", wrap_style, rows = 1 + seq_len(nrow(dictionary)), cols = 3, gridExpand = TRUE)
setColWidths(wb, "Columns", cols = 1:3, widths = c(28, 16, 100))
freezePane(wb, "Columns", firstRow = TRUE)

# One sheet per table
for (sh in tables$sheet) {
  dat <- data_list[[sh]]
  addWorksheet(wb, sh)
  writeData(wb, sh, dat, headerStyle = head_style)
  num_cols <- which(sapply(dat, function(x) is.double(x)))
  if (length(num_cols)) addStyle(wb, sh, num_style, rows = 1 + seq_len(nrow(dat)), cols = num_cols, gridExpand = TRUE)
  setColWidths(wb, sh, cols = seq_along(dat), widths = pmin(60, pmax(10, nchar(names(dat)) + 2,
                                                                     sapply(dat, function(x) max(nchar(as.character(head(x, 200))), na.rm = TRUE)) + 2)))
  freezePane(wb, sh, firstRow = TRUE)
}

saveWorkbook(wb, out_file, overwrite = TRUE)
message("Saved ", out_file, " (", nrow(tables), " tables, ", nrow(dictionary), " column descriptions)")

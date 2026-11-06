# Predict drug response of the PDX tumours with the cell line models, for every feature set and
# correction. PDX is only a test set: every choice below (which model, which ones go in the
# ensemble) is made from the cell line results, never from PDX.
#   - Single best model per drug: of the best model of each dataset (04), the one that did
#     best on the other cell line datasets (mean R-squared in 05). That's the model that
#     generalizes best across labs, which is what PDX asks for.
#   - Ensemble: average of the best models of all datasets, leaving out models with
#     CV R-squared <= 0 (they learned nothing useful).
#
# How the PDX predictions are judged:
#   - R-squared, Pearson and Spearman against the course PDX response. That value is a rescaled
#     best % change in tumour volume, not a dose-response AUC like in cell lines, so it sits on
#     another scale: R-squared is strongly negative for that reason alone, and the correlations
#     are the fair measures.
#   - mRECIST (Xeva, 01): tumours that kept growing under the drug (PD, progressive disease) vs
#     the rest. ROC AUC of the predicted AUC: 0.5 = chance, above 0.5 = higher predicted AUC
#     (more resistant) for the tumours that kept growing.
# Erlotinib (25) and Trametinib (37) have few PDX samples, so their numbers are noisy.
#
# Run from the project root, after 05.

library(data.table)
source("R/feature_filter.R")
source("R/model_utils.R")

drugs <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")
cv    <- fread("results/tables/cv_metrics.csv")[best == TRUE]
cross <- fread("results/tables/cross_dataset_metrics.csv")
xeva  <- readRDS("data/processed/PDX_xeva_response.rds")

# Candidate models and how the single best and the ensemble are chosen
gen  <- cross[, .(cross_r2 = mean(r2)), by = .(feature_set, correction, drug, dataset = train)]
cand <- merge(cv[, .(feature_set, correction, drug, dataset, method, cv_r2 = r2)], gen,
              by = c("feature_set", "correction", "drug", "dataset"))
cand[, single := cross_r2 == max(cross_r2), by = .(feature_set, correction, drug)]
cand[, in_ensemble := cv_r2 > 0]
cand[, in_ensemble := if (any(in_ensemble)) in_ensemble else TRUE, by = .(feature_set, correction, drug)]
fwrite(cand, "results/tables/pdx_model_choice.csv")

preds <- list(); metrics <- list()
for (cfg in split(cand, by = c("feature_set", "correction"))) {
  fs <- cfg$feature_set[1]; cr <- cfg$correction[1]
  pdx <- load_data("PDX_all", fs, cr)

  for (d in drugs) {
    mods <- cfg[drug == d]
    resp <- merge(get_response(pdx, d), xeva[drug == d, .(sample_id, mRECIST, best.average.response)],
                  by = "sample_id", all.x = TRUE)
    p <- sapply(seq_len(nrow(mods)), function(i)
      predict_saved(readRDS(model_path(d, mods$dataset[i], mods$method[i], fs, cr)), pdx, resp$sample_id))
    p <- matrix(p, nrow = nrow(resp), dimnames = list(NULL, mods$dataset))
    out <- list(single = p[, which(mods$single)[1]], ensemble = rowMeans(p[, mods$in_ensemble, drop = FALSE]))
    used <- list(single = mods[single == TRUE, paste0(dataset, " (", method, ")")][1],
                 ensemble = paste(mods[in_ensemble == TRUE, paste0(dataset, " (", method, ")")], collapse = " + "))

    preds[[length(preds) + 1]] <- data.table(feature_set = fs, correction = cr, drug = d, sample_id = resp$sample_id,
                                             obs = resp$value, mRECIST = resp$mRECIST,
                                             single = out$single, ensemble = out$ensemble, p)
    for (m in names(out)) {
      metrics[[length(metrics) + 1]] <- data.table(
        feature_set = fs, correction = cr, drug = d, model = m, models = used[[m]],
        reg_metrics(resp$value, out[[m]]),
        spearman_best_avg_response = cor(out[[m]], resp$best.average.response, method = "spearman", use = "complete.obs"),
        n_mrecist = sum(!is.na(resp$mRECIST)), n_pd = sum(resp$mRECIST == "PD", na.rm = TRUE),
        roc_auc_pd = roc_auc(out[[m]], resp$mRECIST == "PD"))
    }
  }
}
fwrite(rbindlist(preds, fill = TRUE), "results/tables/pdx_predictions.csv")
metrics <- rbindlist(metrics)
fwrite(metrics, "results/tables/pdx_metrics.csv")

print(metrics[, .(r2 = round(mean(r2), 2), pearson = round(mean(pearson), 3), spearman = round(mean(spearman), 3),
                  roc_auc_pd = round(mean(roc_auc_pd), 3)), by = .(feature_set, correction, model)])

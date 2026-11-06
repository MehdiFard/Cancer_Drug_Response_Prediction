# Train one model per drug and cell line dataset with 10-fold CV, for each combination of
#   feature set: multiomics (expression + CNV + mutations) or panel (CNV + mutations of 882 cancer genes)
#   correction:  uncorrected or corrected (batch-corrected in 02)
# with the 3 model types kept in 03. The best of the three per drug and dataset (highest CV
# R-squared) is "the best model" used in 05 and 06.
#
# Same folds (seed 42) and same features in every combination: the feature filter is fitted on
# the uncorrected training samples. (Correction leaves within-dataset variance unchanged, so
# this is the same choice either way.) PDX is not touched here.
#
# Run from the project root, after 02 and 03. Takes about an hour.

library(data.table)
library(caret)
source("R/feature_filter.R")
source("R/model_utils.R")

datasets     <- c("CCLE", "GDSC_broad", "GDSC_sanger")
drugs        <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")
feature_sets <- c("multiomics", "panel")
corrections  <- c("uncorrected", "corrected")
methods      <- fread("results/tables/model_comparison_summary.csv")[kept == TRUE, method]
min_n        <- 20   # Ribociclib has 3 GDSC_sanger samples
n_threads    <- max(1, parallel::detectCores() - 2)
doParallel::registerDoParallel(cores = n_threads)

metrics <- list(); oof <- list()
for (fs in feature_sets) {
  for (ds in datasets) {
    raw <- load_data(ds, fs)
    for (cr in corrections) {
      data <- if (cr == "uncorrected") raw else load_data(ds, fs, cr)
      dir.create(file.path("models", paste0(fs, "_", cr)), recursive = TRUE, showWarnings = FALSE)

      for (drug in drugs) {
        resp <- get_response(raw, drug)
        if (nrow(resp) < min_n) next
        filter <- fit_filter_for(subset_samples(raw, resp$sample_id), fs)
        x <- apply_feature_filter(subset_samples(data, resp$sample_id), filter)
        y <- resp$value

        set.seed(42)
        folds <- createFolds(y, k = 10, returnTrain = TRUE)
        for (m in methods) {
          message(sprintf("%s / %s / %s / %s / %s: %d samples, %d features", fs, cr, drug, ds, m, nrow(x), ncol(x)))
          ctrl <- trainControl(method = "cv", index = folds, savePredictions = "final",
                               allowParallel = !(m %in% threaded_methods))
          set.seed(42)
          fit <- train_model(x, y, m, ctrl, n_threads)

          # Out-of-fold predictions: every sample predicted by a model that never saw it
          pred <- as.data.table(fit$pred)[order(rowIndex)]
          oof[[length(oof) + 1]] <- data.table(feature_set = fs, correction = cr, drug = drug, dataset = ds, method = m,
                                               sample_id = resp$sample_id[pred$rowIndex], obs = pred$obs, pred = pred$pred)
          metrics[[length(metrics) + 1]] <- data.table(
            feature_set = fs, correction = cr, drug = drug, dataset = ds, method = m,
            reg_metrics(pred$obs, pred$pred), n_features = ncol(x),
            tune = paste(names(fit$bestTune), signif(unlist(fit$bestTune), 3), sep = "=", collapse = ", "))

          saveRDS(list(fit = slim_fit(fit), filter = filter, drug = drug, dataset = ds, method = m,
                       feature_set = fs, correction = cr, train_ids = resp$sample_id),
                  model_path(drug, ds, m, fs, cr))
        }
      }
      rm(data); gc()
    }
  }
}
metrics <- rbindlist(metrics)
metrics[, best := r2 == max(r2), by = .(feature_set, correction, drug, dataset)]
fwrite(metrics, "results/tables/cv_metrics.csv")
fwrite(rbindlist(oof), "results/tables/cv_predictions.csv")

print(metrics[best == TRUE, .(feature_set, correction, drug, dataset, method, n, r2 = round(r2, 3))])

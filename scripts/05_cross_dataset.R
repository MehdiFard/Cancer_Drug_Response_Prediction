# Do the models work on data from another lab? The best model of each cell line dataset is
# applied to the other two, for each drug: up to six R-squared values per drug
# (3 training x 2 test datasets), for every feature set and correction.
#
# CCLE and GDSC_broad share most of their cell lines, so those two directions
# aren't really independent tests. n_shared counts the test samples the model already saw
# in training. Tests against GDSC_sanger are the clean ones.
#
# Run from the project root, after 04.

library(data.table)
source("R/feature_filter.R")
source("R/model_utils.R")

datasets <- c("CCLE", "GDSC_broad", "GDSC_sanger")
drugs    <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")
best     <- fread("results/tables/cv_metrics.csv")[best == TRUE]

res <- list()
for (cfg in split(best, by = c("feature_set", "correction"))) {
  fs <- cfg$feature_set[1]; cr <- cfg$correction[1]
  data <- lapply(setNames(datasets, datasets), load_data, feature_set = fs, correction = cr)

  for (i in seq_len(nrow(cfg))) {
    b     <- cfg[i]
    model <- readRDS(model_path(b$drug, b$dataset, b$method, fs, cr))
    for (test_ds in setdiff(datasets, b$dataset)) {
      resp <- get_response(data[[test_ds]], b$drug)
      if (nrow(resp) < 10) next   # R-squared on 3 Ribociclib values means nothing
      pred <- predict_saved(model, data[[test_ds]], resp$sample_id)
      res[[length(res) + 1]] <- data.table(
        feature_set = fs, correction = cr, drug = b$drug, train = b$dataset, test = test_ds, method = b$method,
        reg_metrics(resp$value, pred), n_shared = sum(resp$sample_id %in% model$train_ids))
    }
  }
  rm(data); gc()
}
res <- rbindlist(res)
fwrite(res, "results/tables/cross_dataset_metrics.csv")

# Mean and median per drug, and by training / test dataset (is one dataset giving poor models?)
summ <- res[, .(n_pairs = .N, mean_r2 = mean(r2), median_r2 = median(r2), mean_pearson = mean(pearson, na.rm = TRUE)),
            by = .(feature_set, correction, drug)][order(feature_set, correction, match(drug, drugs))]
fwrite(summ, "results/tables/cross_dataset_summary.csv")
by_ds <- rbind(res[, .(role = "train", mean_r2 = mean(r2), median_r2 = median(r2)), by = .(feature_set, correction, dataset = train)],
               res[, .(role = "test", mean_r2 = mean(r2), median_r2 = median(r2)), by = .(feature_set, correction, dataset = test)])
fwrite(by_ds, "results/tables/cross_dataset_by_dataset.csv")

print(res[, .(mean_r2 = round(mean(r2), 3)), by = .(feature_set, correction)])
print(by_ds[, .(feature_set, correction, role, dataset, mean_r2 = round(mean_r2, 3))])

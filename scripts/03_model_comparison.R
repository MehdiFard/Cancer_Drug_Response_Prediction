# Which model types work best? Six types compared with 10-fold CV on the same folds, for every
# drug and cell line dataset (multi-omics, no batch correction): elastic net (glmnet), random
# forest (ranger), kNN, SVM with radial kernel, PLS and gradient boosted trees (xgboost).
# The 3 with the highest mean CV R-squared are kept and used in 04 (figure in 07).
#
# The training set is the samples with a response for that drug. The feature filter is fitted
# on that training set only (R/feature_filter.R). PDX is not touched here.
#
# Run from the project root. Takes about an hour.

library(data.table)
library(caret)
source("R/feature_filter.R")
source("R/model_utils.R")

datasets  <- c("CCLE", "GDSC_broad", "GDSC_sanger")
drugs     <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")
methods   <- c("glmnet", "ranger", "knn", "svmRadial", "pls", "xgboost")
n_keep    <- 3
min_n     <- 20   # Ribociclib has 3 GDSC_sanger samples, too few to train on
n_threads <- max(1, parallel::detectCores() - 2)

# Most methods run their CV folds in parallel; random forest and xgboost use their own threads
doParallel::registerDoParallel(cores = n_threads)

metrics <- list()
for (ds in datasets) {
  data <- load_data(ds, "multiomics")
  for (drug in drugs) {
    resp <- get_response(data, drug)
    if (nrow(resp) < min_n) next
    train <- subset_samples(data, resp$sample_id)
    x <- apply_feature_filter(train, fit_filter_for(train, "multiomics"))
    y <- resp$value

    # Same folds for every model type, so they're compared on equal terms
    set.seed(42)
    folds <- createFolds(y, k = 10, returnTrain = TRUE)
    for (m in methods) {
      message(sprintf("%s / %s / %s: %d samples, %d features", drug, ds, m, nrow(x), ncol(x)))
      ctrl <- trainControl(method = "cv", index = folds, savePredictions = "final",
                           allowParallel = !(m %in% threaded_methods))
      set.seed(42)
      fit  <- train_model(x, y, m, ctrl, n_threads)
      pred <- as.data.table(fit$pred)
      metrics[[length(metrics) + 1]] <- data.table(drug = drug, dataset = ds, method = m,
                                                   reg_metrics(pred$obs, pred$pred))
    }
  }
  rm(data); gc()
}
metrics <- rbindlist(metrics)
metrics[, rank := frank(-r2), by = .(drug, dataset)]

# Mean CV R-squared over all drug x dataset combinations (and mean rank, as a check)
summary <- metrics[, .(mean_r2 = mean(r2), median_r2 = median(r2), mean_rank = mean(rank),
                       n_best = sum(rank == 1)), by = method][order(-mean_r2)]
summary[, kept := seq_len(.N) <= n_keep]
fwrite(metrics, "results/tables/model_comparison_cv.csv")
fwrite(summary, "results/tables/model_comparison_summary.csv")

print(summary)

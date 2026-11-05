# Task 1.1: train a model for each drug and cell line dataset with 10-fold CV.
# Six model types are compared on the same folds: elastic net, random forest, kNN,
# SVM (radial), PLS and xgboost. The 3 with the highest mean CV R-squared over all
# drugs and datasets are kept; the best of those per drug and dataset is
# "the best model" used in the next steps.
#
# For each model the training set is the samples with a response for that drug.
# The feature filter is fitted on that training set only (R/feature_filter.R).
# PDX is not touched here.
#
# Run from the project root. Takes about an hour.

library(data.table)
library(caret)
library(ggplot2)
source("R/feature_filter.R")
source("R/model_utils.R")

datasets <- c("CCLE", "GDSC_broad", "GDSC_sanger")
drugs    <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")
methods  <- c("glmnet", "ranger", "knn", "svmRadial", "pls", "xgboost")
n_keep   <- 3
min_n    <- 20   # Pitfall 3: Ribociclib has 3 GDSC_sanger samples, too few to train on
n_threads <- max(1, parallel::detectCores() - 2)

dir.create("models/task1", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Most methods run their CV folds in parallel; random forest and xgboost use their own threads
doParallel::registerDoParallel(cores = n_threads)

metrics <- list()
oof     <- list()

for (ds in datasets) {
  data <- readRDS(file.path("data/processed", paste0(ds, "_clean_multiomics.rds")))

  for (drug in drugs) {
    resp <- get_response(data, drug)
    if (nrow(resp) < min_n) {
      message("Skipping ", drug, " / ", ds, ": only ", nrow(resp), " samples")
      next
    }

    train  <- subset_samples(data, resp$sample_id)
    filter <- fit_feature_filter(train)
    x <- apply_feature_filter(train, filter)
    y <- resp$value

    # Same folds for every model type, so they're compared on equal terms
    set.seed(42)
    folds <- createFolds(y, k = 10, returnTrain = TRUE)

    for (m in methods) {
      message(sprintf("%s / %s / %s: %d samples, %d features", drug, ds, m, nrow(x), ncol(x)))
      ctrl <- trainControl(method = "cv", index = folds, savePredictions = "final",
                           allowParallel = !(m %in% threaded_methods))
      set.seed(42)
      fit <- train_model(x, y, m, ctrl, n_threads)

      # Out-of-fold predictions: every sample predicted by a model that never saw it
      pred <- as.data.table(fit$pred)[order(rowIndex)]
      pred[, `:=`(sample_id = resp$sample_id[rowIndex], drug = drug, dataset = ds, method = m)]
      oof[[length(oof) + 1]] <- pred[, .(drug, dataset, method, sample_id, obs, pred)]

      met <- reg_metrics(pred$obs, pred$pred)
      metrics[[length(metrics) + 1]] <- data.table(drug = drug, dataset = ds, method = m, met,
                                                   tune = paste(names(fit$bestTune), signif(unlist(fit$bestTune), 3),
                                                                sep = "=", collapse = ", "))

      # Drop the parts that aren't needed for predicting (see R/model_utils.R)
      fit <- slim_fit(fit)
      saveRDS(list(fit = fit, filter = filter, drug = drug, dataset = ds, method = m,
                   train_ids = resp$sample_id),
              model_path(drug, ds, m))
    }
  }
  rm(data); gc()
}

metrics <- rbindlist(metrics)
oof     <- rbindlist(oof)

# ---- Keep the methods that work best overall ----
# Mean CV R-squared over all drug x dataset combinations (and mean rank, as a check)
metrics[, rank := frank(-r2), by = .(drug, dataset)]
method_summary <- metrics[, .(mean_r2 = mean(r2), median_r2 = median(r2), mean_rank = mean(rank),
                              n_best = sum(rank == 1)), by = method][order(-mean_r2)]
keep <- method_summary$method[seq_len(n_keep)]
method_summary[, kept := method %in% keep]
fwrite(method_summary, "results/tables/task1_cv_method_summary.csv")

# Models of the dropped methods aren't needed anymore
for (m in setdiff(methods, keep)) {
  unlink(metrics[method == m, model_path(drug, dataset, method)])
}

# Best model per drug and dataset = highest CV R-squared among the kept methods
metrics[, kept := method %in% keep]
metrics[, best := kept & r2 == max(r2[kept]), by = .(drug, dataset)]

fwrite(metrics, "results/tables/task1_cv_metrics.csv")
fwrite(oof, "results/tables/task1_cv_predictions.csv")

# ---- CV R-squared table and heatmap of the best models ----
# Row and column order as asked in the task
best <- metrics[best == TRUE]
best[, drug := factor(drug, levels = rev(drugs))]
best[, dataset := factor(dataset, levels = datasets)]

wide <- dcast(best, drug ~ dataset, value.var = "r2")
wide <- wide[match(drugs, drug)]
fwrite(wide, "results/tables/task1_cv_r2_best.csv")

# Fill every drug x dataset cell, so skipped ones show up as "not trained"
grid <- CJ(drug = factor(drugs, levels = rev(drugs)), dataset = factor(datasets, levels = datasets))
hm <- merge(grid, best[, .(drug, dataset, r2, method)], all.x = TRUE)
hm[, label := ifelse(is.na(r2), "not trained\n(n < 20)", sprintf("%.2f\n%s", r2, method))]

lim <- max(abs(hm$r2), na.rm = TRUE)
p <- ggplot(hm, aes(dataset, drug, fill = r2)) +
  geom_tile(colour = "#fcfcfb", linewidth = 1.5) +
  geom_text(aes(label = label), size = 3.6, colour = "#0b0b0b", lineheight = 0.9) +
  scale_fill_gradient2(low = "#e34948", mid = "#f0efec", high = "#2a78d6", midpoint = 0,
                       limits = c(-lim, lim), na.value = "#ffffff", name = "CV R-squared") +
  scale_x_discrete(position = "top") +
  labs(title = "Task 1.1: 10-fold CV R-squared of the best model per drug and dataset",
       subtitle = "R-squared = 1 - SSE/SST on out-of-fold predictions. Model type under the value.",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        plot.title.position = "plot", plot.background = element_rect(fill = "#fcfcfb", colour = NA))
ggsave("results/figures/task1_cv_r2_heatmap.png", p, width = 7.5, height = 6, dpi = 150)

# ---- Method comparison plot ----
cmp <- copy(metrics)
cmp[, method := factor(method, levels = rev(method_summary$method))]
p <- ggplot(cmp, aes(r2, method, colour = kept)) +
  geom_vline(xintercept = 0, colour = "#8a8984", linewidth = 0.4) +
  geom_point(size = 2.2, alpha = 0.55, position = position_jitter(width = 0, height = 0.12, seed = 1)) +
  geom_point(data = method_summary[, .(method = factor(method, levels = levels(cmp$method)), r2 = mean_r2, kept)],
             shape = 23, size = 4.5, fill = "#ffffff", stroke = 1.3) +
  scale_colour_manual(values = c(`TRUE` = "#2a78d6", `FALSE` = "#8a8984"),
                      labels = c(`TRUE` = "kept", `FALSE` = "dropped"), name = NULL) +
  labs(title = "Task 1.1: CV R-squared of each model type",
       subtitle = "One dot per drug x dataset, diamond = mean. The top 3 by mean are kept.",
       x = "CV R-squared", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e5e0", linewidth = 0.3), legend.position = "top",
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))
ggsave("results/figures/task1_cv_method_comparison.png", p, width = 8, height = 5, dpi = 150)

print(method_summary)

print(metrics[, .(drug, dataset, method, n, r2 = round(r2, 3), pearson = round(pearson, 3), best)])

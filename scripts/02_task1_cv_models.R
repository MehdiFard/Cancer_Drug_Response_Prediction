# Task 1.1: train a model for each drug and cell line dataset with 10-fold CV.
# Two model types (elastic net and random forest), the better one by CV R-squared
# is "the best model" for that drug and dataset.
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
methods  <- c("glmnet", "ranger")
min_n    <- 20   # Pitfall 3: Ribociclib has 3 GDSC_sanger samples, too few to train on
n_threads <- max(1, parallel::detectCores() - 2)

dir.create("models/task1", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# glmnet folds run in parallel; random forest uses its own threads instead
doParallel::registerDoParallel(cores = n_threads)

metrics <- list()
oof     <- list()

for (ds in datasets) {
  data <- readRDS(file.path("data/processed", paste0(ds, "_clean.rds")))

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

    # Same folds for both model types, so they're compared on equal terms
    set.seed(42)
    folds <- createFolds(y, k = 10, returnTrain = TRUE)

    for (m in methods) {
      message(sprintf("%s / %s / %s: %d samples, %d features", drug, ds, m, nrow(x), ncol(x)))
      ctrl <- trainControl(method = "cv", index = folds, savePredictions = "final",
                           allowParallel = (m == "glmnet"))
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

      # The training matrix is big and not needed for predicting
      # (glmnet also keeps a copy inside its call)
      fit$trainingData <- NULL
      fit$finalModel$call <- NULL
      saveRDS(list(fit = fit, filter = filter, drug = drug, dataset = ds, method = m,
                   train_ids = resp$sample_id),
              model_path(drug, ds, m))
    }
  }
  rm(data); gc()
}

metrics <- rbindlist(metrics)
oof     <- rbindlist(oof)

# Best model per drug and dataset = higher CV R-squared
metrics[, best := r2 == max(r2), by = .(drug, dataset)]

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

print(metrics[, .(drug, dataset, method, n, r2 = round(r2, 3), pearson = round(pearson, 3), best)])

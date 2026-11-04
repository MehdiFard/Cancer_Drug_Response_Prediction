# Task 1.2: apply the best model of each cell line dataset to the other two.
# Gives up to six R-squared values per drug (3 training sets x 2 test sets).
#
# Pitfall 5: CCLE and GDSC_broad share most of their cell lines, so those two
# directions aren't really independent tests. The n_shared column counts the test
# samples the model already saw in training. Tests against GDSC_sanger are the clean ones.
#
# Run from the project root, after 02_task1_cv_models.R.

library(data.table)
library(ggplot2)
source("R/feature_filter.R")
source("R/model_utils.R")

datasets <- c("CCLE", "GDSC_broad", "GDSC_sanger")
drugs    <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")
pal      <- c(CCLE = "#2a78d6", GDSC_broad = "#eb6834", GDSC_sanger = "#1baf7a")

cv   <- fread("results/tables/task1_cv_metrics.csv")
best <- cv[best == TRUE]
data <- lapply(setNames(datasets, datasets),
               function(ds) readRDS(file.path("data/processed", paste0(ds, "_clean.rds"))))

res <- list()
for (i in seq_len(nrow(best))) {
  b     <- best[i]
  model <- readRDS(model_path(b$drug, b$dataset, b$method))

  for (test_ds in setdiff(datasets, b$dataset)) {
    resp <- get_response(data[[test_ds]], b$drug)
    # Pitfall 3: Ribociclib has 3 GDSC_sanger samples, R-squared on 3 values means nothing
    if (nrow(resp) < 10) next
    pred <- predict_saved(model, data[[test_ds]], resp$sample_id)
    res[[length(res) + 1]] <- data.table(
      drug = b$drug, train = b$dataset, test = test_ds, method = b$method,
      reg_metrics(resp$value, pred),
      n_shared = sum(resp$sample_id %in% model$train_ids))
  }
}
res <- rbindlist(res)
fwrite(res, "results/tables/task1_cross_dataset_metrics.csv")

# Mean and median R-squared per drug
summ <- res[, .(n_pairs = .N, mean_r2 = mean(r2), median_r2 = median(r2),
                mean_pearson = mean(pearson, na.rm = TRUE)), by = drug]
summ <- summ[match(drugs, drug)]
fwrite(summ, "results/tables/task1_cross_dataset_summary.csv")

# Is one dataset giving poor models? Average R-squared by training set and by test set.
by_train <- res[, .(mean_r2 = mean(r2), median_r2 = median(r2)), by = .(dataset = train)][, role := "train"]
by_test  <- res[, .(mean_r2 = mean(r2), median_r2 = median(r2)), by = .(dataset = test)][, role := "test"]
fwrite(rbind(by_train, by_test), "results/tables/task1_cross_dataset_by_dataset.csv")

# ---- Boxplot of R-squared per drug ----
res[, drug := factor(drug, levels = drugs)]
p <- ggplot(res, aes(drug, r2)) +
  geom_hline(yintercept = 0, colour = "#8a8984", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.4, colour = "#52514e", fill = "#f0efec") +
  geom_point(aes(colour = train, shape = test), size = 2.8, stroke = 0.8,
             position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  scale_colour_manual(values = pal, name = "Trained on") +
  scale_shape_manual(values = c(CCLE = 16, GDSC_broad = 17, GDSC_sanger = 15), name = "Tested on") +
  labs(title = "Task 1.2: models applied to the other cell line datasets",
       subtitle = "R-squared = 1 - SSE/SST. Below 0 = worse than predicting the test set's mean.",
       x = NULL, y = "R-squared") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e5e0", linewidth = 0.3),
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))
ggsave("results/figures/task1_cross_dataset_r2_boxplot.png", p, width = 10, height = 5.5, dpi = 150)

print(res[, .(drug, train, test, method, n, n_shared, r2 = round(r2, 3), pearson = round(pearson, 3))])
print(summ)
print(rbind(by_train, by_test))

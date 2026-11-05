# Task 1.3: predict PDX drug response with the single best model and with an ensemble.
#
# PDX is only a test set: every choice below (which model, which ones go in the
# ensemble) is made from the cell line results, never from PDX.
#   - Single best model per drug: of the best model of each dataset (Task 1.1), the one
#     that did best on the other cell line datasets (mean R-squared in Task 1.2).
#     That's the model that generalizes best across labs, which is what PDX asks for.
#   - Ensemble: average of the best models of all datasets, leaving out models with
#     CV R-squared <= 0 (they learned nothing useful).
#
# PDX AUC comes from tumor volume, not cell viability, so the values sit on another
# scale (course text + our figures). R-squared (1 - SSE/SST) will be strongly negative
# for that reason alone; Pearson and Spearman are the fair measures here.
# Pitfall 2: Erlotinib (25) and Trametinib (37) have few PDX samples, so their
# correlations are noisy. Pitfall 8: PDX CNV is almost all 0.
#
# Run from the project root, after 02 and 03.

library(data.table)
library(ggplot2)
source("R/feature_filter.R")
source("R/model_utils.R")

drugs <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")

cv    <- fread("results/tables/task1_cv_metrics.csv")[best == TRUE]
cross <- fread("results/tables/task1_cross_dataset_metrics.csv")
pdx   <- readRDS("data/processed/PDX_all_clean_multiomics.rds")

# Pick the single best model per drug from the cell line results
gen <- cross[, .(cross_r2 = mean(r2)), by = .(drug, dataset = train)]
cand <- merge(cv[, .(drug, dataset, method, cv_r2 = r2)], gen, by = c("drug", "dataset"))
cand[, single := cross_r2 == max(cross_r2), by = drug]
cand[, in_ensemble := cv_r2 > 0]
cand[, in_ensemble := if (any(in_ensemble)) in_ensemble else TRUE, by = drug]  # never empty
fwrite(cand, "results/tables/task1_pdx_model_choice.csv")

preds   <- list()
metrics <- list()
for (d in drugs) {
  resp <- get_response(pdx, d)
  mods <- cand[drug == d]

  # Prediction from every candidate model
  p <- sapply(seq_len(nrow(mods)), function(i) {
    predict_saved(readRDS(model_path(d, mods$dataset[i], mods$method[i])), pdx, resp$sample_id)
  })
  p <- matrix(p, nrow = nrow(resp), dimnames = list(NULL, mods$dataset))

  single   <- p[, which(mods$single)[1]]
  ensemble <- rowMeans(p[, mods$in_ensemble, drop = FALSE])

  preds[[d]] <- data.table(drug = d, sample_id = resp$sample_id, obs = resp$value,
                           single = single, ensemble = ensemble, p)
  metrics[[d]] <- rbind(
    data.table(drug = d, model = "single", models = paste(mods[single == TRUE, paste0(dataset, " (", method, ")")]),
               reg_metrics(resp$value, single)),
    data.table(drug = d, model = "ensemble", models = paste(mods[in_ensemble == TRUE, paste0(dataset, " (", method, ")")], collapse = " + "),
               reg_metrics(resp$value, ensemble)))
}
preds   <- rbindlist(preds, fill = TRUE)
metrics <- rbindlist(metrics)
fwrite(preds, "results/tables/task1_pdx_predictions.csv")
fwrite(metrics, "results/tables/task1_pdx_metrics.csv")

# ---- Plot: R-squared, Pearson and Spearman per drug, single vs ensemble ----
long <- melt(metrics, id.vars = c("drug", "model", "n"), measure.vars = c("r2", "pearson", "spearman"),
             variable.name = "metric")
long[, metric := factor(metric, levels = c("r2", "pearson", "spearman"),
                        labels = c("R-squared (1 - SSE/SST)", "Pearson r", "Spearman rho"))]
n_lab <- metrics[model == "single", setNames(sprintf("%s\nn = %d", drug, n), drug)]
long[, drug := factor(drug, levels = drugs, labels = n_lab[drugs])]
long[, model := factor(model, levels = c("single", "ensemble"), labels = c("Single best model", "Ensemble (average)"))]

p <- ggplot(long, aes(drug, value, fill = model)) +
  geom_hline(yintercept = 0, colour = "#8a8984", linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "#fcfcfb", linewidth = 0.5) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Single best model" = "#2a78d6", "Ensemble (average)" = "#eb6834"), name = NULL) +
  labs(title = "Task 1.3: predicting PDX drug response from models trained on cell lines",
       subtitle = "PDX AUC is on another scale than cell line AUC, so R-squared is negative by construction; correlations are the fair test.\nPositive correlation = models rank PDX samples in the right order.",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e5e0", linewidth = 0.3),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "top",
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))
ggsave("results/figures/task1_pdx_metrics.png", p, width = 10, height = 8.5, dpi = 150)

print(cand)
print(metrics[, .(drug, model, n, r2 = round(r2, 2), pearson = round(pearson, 3), spearman = round(spearman, 3), models)])

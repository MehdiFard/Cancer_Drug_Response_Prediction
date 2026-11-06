# Which features drive the predictions? Read from the elastic net models (multi-omics,
# batch-corrected), one per drug and cell line dataset, trained in 04. Elastic net is the one of the
# three kept model types with a weight per feature that can be read directly; SVM has none.
#
# Weights are standardized: coefficient x SD of the feature in the training samples, i.e. the
# change in predicted AUC for a 1 SD change in the feature. Positive = higher value, higher
# predicted AUC = more resistant. Negative = more sensitive.
# Per drug they are averaged over the datasets (0 where a dataset's model didn't pick the
# feature), and n_datasets counts how many models picked it: a feature picked in every dataset
# is a more reliable marker than one picked once.
#
# Saves results/tables/top_features.csv (top 20 per drug). Run from the project root, after 04.

library(data.table)
source("R/feature_filter.R")
source("R/model_utils.R")

drugs    <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")
datasets <- c("CCLE", "GDSC_broad", "GDSC_sanger")
n_top    <- 20

w <- list()
for (ds in datasets) {
  data <- load_data(ds, "multiomics", "corrected")
  for (d in drugs) {
    path <- model_path(d, ds, "glmnet", "multiomics", "corrected")
    if (!file.exists(path)) next   # Ribociclib in GDSC_sanger: not trained
    m    <- readRDS(path)
    x    <- apply_feature_filter(subset_samples(data, m$train_ids), m$filter)
    beta <- as.matrix(glmnet::coef.glmnet(m$fit$finalModel, s = m$fit$bestTune$lambda))[-1, 1]   # drop intercept
    w[[length(w) + 1]] <- data.table(drug = d, dataset = ds, feature = names(beta),
                                     weight = beta * apply(x[, names(beta), drop = FALSE], 2, sd))
  }
  rm(data); gc()
}
w <- rbindlist(w)
w[, n_models := uniqueN(dataset), by = drug]

top <- w[, .(mean_weight = sum(weight) / n_models[1], n_datasets = sum(weight != 0), n_models = n_models[1]),
         by = .(drug, feature)]
top <- top[n_datasets > 0][order(drug, -abs(mean_weight))][, head(.SD, n_top), by = drug]

# Feature names are <omic>_<Ensembl ID>; add gene symbols
top[, `:=`(omic = sub("_.*", "", feature), ensembl = sub("^[a-z]+_", "", feature))]
sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(top$ensembl),
                                              column = "SYMBOL", keytype = "ENSEMBL"))
top[, gene := ifelse(is.na(sym[ensembl]), ensembl, sym[ensembl])]
top[, rank := seq_len(.N), by = drug]
setcolorder(top, c("drug", "rank", "gene", "omic", "mean_weight", "n_datasets", "n_models", "ensembl", "feature"))
fwrite(top[order(match(drug, drugs), rank)], "results/tables/top_features.csv")

print(top[rank <= 5, .(drug, rank, gene, omic, mean_weight = signif(mean_weight, 2), n_datasets)])

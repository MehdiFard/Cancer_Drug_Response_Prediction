# Batch correction across the four datasets, for both feature sets.
# The datasets were measured in different labs and on different platforms, and that shows
# up as large systematic differences (e.g. CCLE vs GDSC, although they share ~540 cell lines).
# regressBatches removes each dataset's own mean per feature (R/batch_correction.R).
# This uses the PDX omics only, never the PDX responses.
#
# Every omic is corrected separately. Mutations (0/1) and CNV (-1/0/1) get the same mean
# removal; setting mutation 0s to -1 first would change nothing, since removing the mean
# is the same on 0/1 and -1/1 data up to a factor 2.
# PDX CNV is ~99% 0 (much stricter calls than in cell lines), and no correction can bring back
# information that isn't there.
#
# Saves data/processed/<dataset>_corrected_<feature set>.rds, and PCA coordinates before / after
# (results/tables/batch_effect_pca.csv, drawn in 08).
# Run from the project root. Takes about a minute.

library(data.table)
source("R/batch_correction.R")

datasets <- c("CCLE", "GDSC_broad", "GDSC_sanger", "PDX_all")

for (fs in c("multiomics", "panel")) {
  data  <- lapply(setNames(datasets, datasets),
                  function(ds) readRDS(file.path("data/processed", paste0(ds, "_clean_", fs, ".rds"))))
  omics <- setdiff(names(data[[1]]), "drug")
  corrected <- lapply(setNames(omics, omics), function(o) regress_batches(lapply(data, function(x) x[[o]] * 1.0)))
  for (ds in datasets) {
    out <- c(lapply(corrected, `[[`, ds), list(drug = data[[ds]]$drug))
    saveRDS(out, file.path("data/processed", paste0(ds, "_corrected_", fs, ".rds")))
  }
  if (fs == "multiomics") gex <- list(uncorrected = lapply(data, `[[`, "gex"), corrected = corrected$gex)
}

# ---- PCA of expression before and after ----
# The same 2000 genes in both: the most variable over all samples pooled, before correction.
# These include the genes that differ most between datasets, so the PCA shows the batch effect
# and what the correction does to it. Genes are scaled to unit variance.
v   <- apply(do.call(rbind, gex$uncorrected), 2, var)
top <- names(sort(v, decreasing = TRUE))[1:2000]

pcs <- list(); mixing <- list()
for (cr in names(gex)) {
  X  <- do.call(rbind, lapply(gex[[cr]], function(m) m[, top]))
  ds <- rep(datasets, sapply(gex[[cr]], nrow))
  pc <- prcomp(X, center = TRUE, scale. = TRUE, rank. = 10)

  # How mixed are the datasets? Share of each sample's 20 nearest neighbours (top 10 PCs)
  # that come from its own dataset. Fully mixed would give ~sum(p^2), p = dataset shares.
  nn   <- BiocNeighbors::findKNN(pc$x, k = 20)$index
  same <- rowMeans(matrix(ds[nn], nrow = nrow(nn)) == ds)
  mixing[[cr]] <- data.table(correction = cr,
                             same_dataset_neighbours_cell_lines = mean(same[ds != "PDX_all"]),
                             same_dataset_neighbours_pdx = mean(same[ds == "PDX_all"]),
                             if_fully_mixed = sum(prop.table(table(ds))^2))
  pcs[[cr]] <- data.table(correction = cr, sample_id = rownames(X), dataset = ds,
                          PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                          var_PC1 = summary(pc)$importance[2, 1], var_PC2 = summary(pc)$importance[2, 2])
}
mixing <- rbindlist(mixing)
fwrite(mixing, "results/tables/batch_effect_mixing.csv")
fwrite(rbindlist(pcs), "results/tables/batch_effect_pca.csv")   # drawn in 08

print(mixing)

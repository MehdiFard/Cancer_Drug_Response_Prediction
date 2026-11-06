# Feature filtering.
# Which features to keep is learned on the training set only and then applied to the
# test sets, so the test data never influences the choice.
#
# Expression and CNV: top 5000 most variable genes.
# Mutations: genes mutated in at least 5% of samples. 1% keeps ~13,000 genes in the
# multi-omics data, mostly long genes hit by chance; 5% keeps ~1,600-2,800.
# The gene panel has no expression and only curated cancer genes, so there all CNV genes
# are kept and mutations down to 1% (e.g. EGFR mutations, rare in cell lines but key for Erlotinib).

# Optional correlation filter (off by default): set cor_cutoff, e.g. 0.9, to drop one gene
# from each highly correlated pair in expression and CNV. CNV has lots of these,
# since neighbouring genes usually share the same copy number segment.

# Learn which features to keep from one training dataset (a list with cnv, mut and maybe gex)
fit_feature_filter <- function(train, n_gex = 5000, n_cnv = 5000, min_mut_freq = 0.05,
                               cor_cutoff = NULL) {

  # Drop near-zero-variance columns, keep the n with the highest SD,
  # then (if asked) drop highly correlated ones
  top_variable <- function(m, n) {
    nzv  <- colnames(m)[caret::nearZeroVar(m, uniqueCut = 15)]
    m    <- m[, setdiff(colnames(m), nzv), drop = FALSE]
    sds  <- apply(m, 2, sd)
    keep <- names(sort(sds, decreasing = TRUE))[seq_len(min(n, length(sds)))]

    if (!is.null(cor_cutoff)) {
      too_similar <- caret::findCorrelation(cor(m[, keep]), cutoff = cor_cutoff, names = TRUE)
      keep <- setdiff(keep, too_similar)
    }
    keep
  }

  filter <- list()
  if (!is.null(train$gex)) filter$gex <- top_variable(train$gex, n_gex)
  filter$cnv <- top_variable(train$cnv, n_cnv)
  # nzv would remove almost every mutation, so we use mutation frequency instead
  filter$mut <- colnames(train$mut)[colMeans(train$mut) >= min_mut_freq]
  filter
}

# Settings per feature set
fit_filter_for <- function(train, feature_set) {
  switch(feature_set,
         multiomics = fit_feature_filter(train),
         panel      = fit_feature_filter(train, min_mut_freq = 0.01))
}

# Apply a fitted filter to any dataset and glue the omics into one matrix.
# Column names get a prefix, so the same gene can appear once per omic.
apply_feature_filter <- function(data, filter) {
  x <- do.call(cbind, lapply(names(filter), function(o) data[[o]][, filter[[o]], drop = FALSE]))
  colnames(x) <- unlist(lapply(names(filter), function(o) paste0(o, "_", filter[[o]])))
  x
}

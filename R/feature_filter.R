# Feature filtering for the multi-omics matrix.
# Which features to keep is learned on the training set only and then applied to the
# test sets, so the test data never influences the choice (the course is strict about this).
#
# Expression and CNV: top 5000 most variable genes, as in 2021.
# Mutations: genes mutated in at least 5% of samples. 1% (≈ the 2021 SD > 0.1 rule) keeps
# ~13,000 genes at gene level, mostly long genes hit by chance; 5% keeps ~1,600-2,800.

# Optional correlation filter (off by default): set cor_cutoff, e.g. 0.9, to drop one gene
# from each highly correlated pair in expression and CNV (as in 2021). CNV has lots of these,
# since neighbouring genes usually share the same copy number segment.

# Learn which features to keep from one training dataset (a list with gex, cnv, mut)
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

  list(
    gex = top_variable(train$gex, n_gex),
    cnv = top_variable(train$cnv, n_cnv),
    # nzv would remove almost every mutation, so we use mutation frequency instead
    mut = colnames(train$mut)[colMeans(train$mut) >= min_mut_freq]
  )
}

# Apply a fitted filter to any dataset and glue the omics into one matrix.
# Column names get a prefix, so the same gene can appear once per omic.
apply_feature_filter <- function(data, filter) {
  x <- cbind(data$gex[, filter$gex, drop = FALSE],
             data$cnv[, filter$cnv, drop = FALSE],
             data$mut[, filter$mut, drop = FALSE])
  colnames(x) <- c(paste0("gex_", filter$gex),
                   paste0("cnv_", filter$cnv),
                   paste0("mut_", filter$mut))
  x
}

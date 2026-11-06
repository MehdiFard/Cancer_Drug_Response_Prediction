# Batch correction across datasets.
# Each dataset (CCLE, GDSC_broad, GDSC_sanger, PDX_all) is one batch. Takes a named list of
# sample x feature matrices, one per dataset, and returns the same list corrected.
#
# regressBatches (batchelor) fits feature ~ batch and keeps the residuals. With only the batch
# in the model that's the same as removing each dataset's own mean per feature (checked:
# identical). So each dataset is corrected with its own samples only: adding PDX doesn't change
# the corrected cell line data, and a new batch of samples can be corrected the same way later.
# Per-dataset z-scoring and fastMNN were tried as well and did no better.
regress_batches <- function(mats) {
  batch <- rep(names(mats), sapply(mats, nrow))
  out <- batchelor::regressBatches(t(do.call(rbind, mats)), batch = batch)
  out <- t(as.matrix(SummarizedExperiment::assay(out, "corrected")))

  # Split back by position, since CCLE and GDSC_broad share sample IDs
  lapply(setNames(names(mats), names(mats)), function(s) {
    m <- out[batch == s, , drop = FALSE]
    dimnames(m) <- dimnames(mats[[s]])
    m
  })
}

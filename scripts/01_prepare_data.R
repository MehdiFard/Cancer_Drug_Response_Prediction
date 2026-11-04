# Task 1: read and prepare the data
# Reads expression, CNV, mutations and drug response for all four datasets,
# fixes the known data problems and saves one clean .rds per dataset in data/processed/.
#
# Only steps that don't depend on the training set happen here. Feature filtering
# (top variable genes etc.) has to be learned on the training data only, so it lives
# in R/feature_filter.R and runs later, when we train.
#
# Run from the project root.

library(data.table)

raw_dir  <- "data/raw"
out_dir  <- "data/processed"
datasets <- c("CCLE", "GDSC_broad", "GDSC_sanger", "PDX_all")
drugs    <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")


# ---- Helpers ----

# gex and cnv are wide (samples x genes), sample IDs sit in the first column.
# GDSC files are actually comma-separated despite the .tsv name, fread sorts that out.
read_wide <- function(path) {
  dt <- fread(path)
  m  <- as.matrix(dt[, -1])
  rownames(m) <- dt[[1]]
  m
}

# One file per drug, stacked into one long table: sample_id, drug, value.
# The other two columns are just labels (and swapped in the PDX files), so we skip them.
read_drugs <- function(dataset) {
  files <- file.path(raw_dir, dataset, "drug_response", paste0(drugs, ".tsv.gz"))
  dt <- rbindlist(lapply(files, fread,
                         select = c("sample_id", "column_name", "value"),
                         colClasses = list(numeric = "value")))
  setnames(dt, "column_name", "drug")
  dt
}

# Mutations come as one row per mutation. We only need sample, gene and protein change.
# Pitfall 7: PDX has its own format, the gene ID is in "Gene" and protein changes
# are written without the "p." (C1277G instead of p.C1277G).
read_mut <- function(dataset) {
  path <- file.path(raw_dir, dataset, "all_intersect_genes", "mut.tsv.gz")
  if (dataset == "PDX_all") {
    dt <- fread(path, select = c("sample_id", "Gene", "protein_change"))
    setnames(dt, "Gene", "gene_id")
  } else {
    dt <- fread(path, select = c("sample_id", "gene_id", "protein_change"))
  }
  dt
}

# Gene-level 0/1 matrix: 1 if the sample has at least one non-silent mutation in that gene.
# Silent = same amino acid before and after (p.A2A). Only CCLE labels them (Variant_Classification),
# but this pattern matches that label almost exactly and works for every dataset.
# Using the full gene list keeps the columns identical across datasets.
binarize_mut <- function(mut, samples, genes) {
  aa_change <- sub("^p\\.", "", mut$protein_change)
  silent    <- grepl("^([A-Z])[0-9]+\\1$", aa_change)
  hits <- unique(mut[!silent & sample_id %in% samples & gene_id %in% genes, .(sample_id, gene_id)])

  m <- matrix(0L, nrow = length(samples), ncol = length(genes), dimnames = list(samples, genes))
  m[cbind(hits$sample_id, hits$gene_id)] <- 1L
  m
}

# Fill missing values from the most similar samples (k nearest neighbours, within one dataset).
# caret's knnImpute also centers and scales every column, so we undo that afterwards
# to keep the original values.
knn_impute <- function(m) {
  pp  <- caret::preProcess(m, method = "knnImpute")
  out <- predict(pp, m)
  sweep(sweep(out, 2, pp$std, "*"), 2, pp$mean, "+")
}


# ---- Read and prepare each dataset ----

dir.create(out_dir, showWarnings = FALSE)

for (ds in datasets) {
  message("Preparing ", ds, " ...")
  omics_dir <- file.path(raw_dir, ds, "all_intersect_genes")

  gex  <- read_wide(file.path(omics_dir, "gex.tsv.gz"))
  cnv  <- read_wide(file.path(omics_dir, "cnv.tsv.gz"))
  mut  <- read_mut(ds)
  drug <- read_drugs(ds)

  # Keep samples that have all three omics, in the same order everywhere
  samples <- sort(Reduce(intersect, list(rownames(gex), rownames(cnv), unique(mut$sample_id))))
  gex <- gex[samples, ]
  cnv <- cnv[samples, ]

  # Pitfall 1: PDX expression is raw FPKM (up to ~325,000). The cell lines are on a
  # log2(x + 1) scale (0 to ~16; the course says log10, but the values don't fit that).
  # Same transform for PDX puts it on the same scale.
  if (ds == "PDX_all") gex <- log2(gex + 1)

  # CNV is already -1 / 0 / 1, nothing to do

  # This data has no missing values, but if some ever show up, impute them.
  # CNV gets rounded back to -1 / 0 / 1 afterwards.
  if (anyNA(gex)) {
    message("  imputing ", sum(is.na(gex)), " missing gex values")
    gex <- knn_impute(gex)
  }
  if (anyNA(cnv)) {
    message("  imputing ", sum(is.na(cnv)), " missing cnv values")
    cnv <- round(knn_impute(cnv))
  }

  mut <- binarize_mut(mut, samples, genes = colnames(gex))

  # Pitfall 2: most PDX drug values are empty, drop those rows.
  # Pitfall 6: some GDSC cell lines were measured twice for Erlotinib and Trametinib.
  # We take the mean. For Trametinib the two values often disagree a lot, see the notes.
  # Pitfall 4: both GDSC folders have the same drug files, but each dataset only keeps
  # its own cell lines here, so that's fine.
  drug <- drug[!is.na(value) & sample_id %in% samples,
               .(value = mean(value)), by = .(sample_id, drug)]

  # Last safety check before saving. Mutations can't be missing after binarizing,
  # and drug responses are never imputed (they're what we predict).
  stopifnot(!anyNA(gex), !anyNA(cnv), !anyNA(mut), !anyNA(drug))

  saveRDS(list(gex = gex, cnv = cnv, mut = mut, drug = drug),
          file.path(out_dir, paste0(ds, "_clean.rds")))
}


# ---- Quick check of what we saved ----

clean <- lapply(setNames(datasets, datasets),
                function(ds) readRDS(file.path(out_dir, paste0(ds, "_clean.rds"))))

# Samples and genes per omic (should be the same within each dataset)
cat("\n== Samples, and genes per omic (genes = expression) ==\n")
sapply(clean, function(x) c(samples = nrow(x$gex), genes = ncol(x$gex),
                            cnv = ncol(x$cnv), mut = ncol(x$mut)))

# Expression should now be on a similar scale everywhere
cat("\n== Expression values, log2(x + 1): min, median, 99th percentile, max ==\n")
sapply(clean, function(x) round(quantile(x$gex, c(0, 0.5, 0.99, 1)), 2))

# Mutated genes per sample (median)
cat("\n== Median number of mutated genes per sample (non-silent) ==\n")
sapply(clean, function(x) median(rowSums(x$mut)))

# Samples with a response, per drug.
# Pitfall 3: Ribociclib has very few cell lines. Pitfall 2: Erlotinib and Trametinib
# have very few PDX samples. Keep both in mind when reading the results.
cat("\n== Samples with a drug response, per drug ==\n")
dcast(rbindlist(lapply(clean, `[[`, "drug"), idcol = "dataset"),
      drug ~ dataset, fun.aggregate = length)

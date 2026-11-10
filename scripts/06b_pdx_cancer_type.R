# Does cancer type explain the PDX results? The course data has no cancer type; Xeva's model.tsv
# has one per PDX sample (breast, melanoma, colorectal, lung, pancreatic, gastric).
#   - Which cancer types each drug was tested on: a drug tested on one cancer type only can't be
#     predicted well just by telling cancer types apart.
#   - Share of the variance explained by cancer type, in the predictions and in the observed PDX
#     response. If it is high in the predictions but low in the response, the models carry
#     cancer type differences from the cell lines that don't hold in PDX.
#   - Spearman over all PDX vs within cancer types (each cancer type's mean rank removed).
# All on ranks, since the PDX response and the predicted AUC sit on different scales.
#
# Saves results/tables/pdx_cancer_type_counts.csv and pdx_cancer_type_stats.csv.
# Run from the project root, after 06.

library(data.table)

xd <- file.path(tempdir(), "xeva")
untar("data/raw/xeva/Xeva_PDXE_tsv.tar.gz", files = "Xeva_PDXE_tsv/model.tsv", exdir = xd)
types <- unique(fread(file.path(xd, "Xeva_PDXE_tsv", "model.tsv"))[, .(sample_id = patient.id, cancer_type = tissue)])
stopifnot(!anyDuplicated(types$sample_id))

pp <- fread("results/tables/pdx_predictions.csv", na.strings = c("", "NA"))[!is.na(obs)]
pp <- merge(pp, types, by = "sample_id", all.x = TRUE)
pp[is.na(cancer_type), cancer_type := "unknown"]   # one Trametinib PDX is not in Xeva

counts <- pp[feature_set == "multiomics" & correction == "corrected",
             .(n = .N, n_pd = sum(mRECIST == "PD", na.rm = TRUE)), by = .(drug, cancer_type)]
fwrite(counts[order(drug, cancer_type)], "results/tables/pdx_cancer_type_counts.csv")

eta2 <- function(y, g) if (uniqueN(g) > 1) summary(lm(y ~ g))$r.squared else NA_real_
pp[, `:=`(rp = frank(ensemble), ro = frank(obs)), by = .(feature_set, correction, drug)]
pp[, n_type := .N, by = .(feature_set, correction, drug, cancer_type)]
stats <- pp[, .(n = .N, n_types = uniqueN(cancer_type[n_type >= 5]),   # cancer types with >= 5 PDX
                var_pred = eta2(rp, cancer_type), var_obs = eta2(ro, cancer_type),
                spearman_all = cor(rp, ro),
                spearman_within = cor(rp - ave(rp, cancer_type), ro - ave(ro, cancer_type))),
            by = .(feature_set, correction, drug)]
fwrite(stats, "results/tables/pdx_cancer_type_stats.csv")

print(dcast(counts, drug ~ cancer_type, value.var = "n", fill = 0))
print(stats[feature_set == "multiomics" & correction == "corrected"], digits = 2)

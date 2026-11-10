# Figures and the overall summary table. Everything is drawn from results/tables/, so the
# figures can be redrawn in seconds without rerunning the models. Captions are in README.md.
#   fig1_data_and_cv      batch effect (PCA), model types, CV R-squared per drug and dataset
#   fig2_cross_dataset    best models applied to the other cell line datasets
#   fig3_pdx              PDX: correlation, mRECIST ROC AUC, predictions by response class, cancer type
#   fig4_top_features     strongest features of the elastic net models per drug
#   summary.csv           one line per feature set and correction
#
# The figures show the batch-corrected results only (the PCA shows both). summary.csv has both,
# for comparison: uncorrected gives the same CV and nearly the same everything else.
#
# Run from the project root, after 06, 06b and 07.

library(data.table)
library(ggplot2)
library(patchwork)

drugs    <- c("ALPELISIB", "BUPARLISIB", "CETUXIMAB", "ERLOTINIB", "RIBOCICLIB", "TRAMETINIB")
datasets <- c("CCLE", "GDSC_broad", "GDSC_sanger")
drug_lab <- setNames(tools::toTitleCase(tolower(drugs)), drugs)
ds_lab   <- c(CCLE = "CCLE", GDSC_broad = "GDSC-Broad", GDSC_sanger = "GDSC-Sanger", PDX_all = "PDX")
fs_lab   <- c(multiomics = "Multi-omics", panel = "Gene panel")
fs_pal   <- setNames(c("#1c5cab", "#c0461b"), fs_lab)
ds_pal   <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#6250d6"), ds_lab)
ink      <- "#1f1f1e"; grey <- "#8a8984"

cv    <- fread("results/tables/cv_metrics.csv")
cross <- fread("results/tables/cross_dataset_metrics.csv")
pdx   <- fread("results/tables/pdx_metrics.csv")

# ---- Summary tables ----
summary <- merge(
  cv[best == TRUE, .(mean_cv_r2 = mean(r2)), by = .(feature_set, correction)],
  cross[, .(mean_cross_dataset_r2 = mean(r2)), by = .(feature_set, correction)], by = c("feature_set", "correction"))
pdx_avg <- pdx[, .(pdx_r2 = mean(r2), pdx_pearson = mean(pearson, na.rm = TRUE), pdx_spearman = mean(spearman, na.rm = TRUE),
                   pdx_roc_auc_pd = mean(roc_auc_pd, na.rm = TRUE)), by = .(feature_set, correction, model)]
pdx_avg <- dcast(pdx_avg, feature_set + correction ~ model, value.var = c("pdx_r2", "pdx_pearson", "pdx_spearman", "pdx_roc_auc_pd"))
summary <- merge(summary, pdx_avg, by = c("feature_set", "correction"))
fwrite(summary, "results/tables/summary.csv")

# Best CV R-squared per drug and dataset, one column per feature set and correction
cv_wide <- dcast(cv[best == TRUE], drug + dataset ~ feature_set + correction, value.var = "r2")
fwrite(cv_wide[order(match(drug, drugs), dataset)], "results/tables/cv_r2_best.csv")

# ---- Style ----
# Journal style: small sans-serif text, axis lines without a box, no titles (captions in README),
# lowercase bold panel letters. Width 180 mm = a two-column figure.
th <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink, size = 7),
        axis.line = element_line(linewidth = 0.3, colour = ink), axis.ticks = element_line(linewidth = 0.3, colour = ink),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7), legend.title = element_text(size = 7),
        legend.background = element_blank(), plot.tag = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(face = "bold", size = 8, hjust = 0.5),
        plot.margin = margin(3, 5, 3, 3))
tags <- plot_annotation(tag_levels = "a")

save_fig <- function(p, name, width, height) {
  ggsave(file.path("results/figures", paste0(name, ".png")), p, width = width, height = height, units = "mm",
         dpi = 300, device = ragg::agg_png, bg = "white")
}

# ---- Figure 1: data and cross-validation ----
# a-b batch effect
pcs <- fread("results/tables/batch_effect_pca.csv")
pcs[, dataset := factor(ds_lab[dataset], levels = ds_lab)]
set.seed(1)

pca_panel <- function(cr) {
  d <- pcs[correction == cr]
  ggplot(d[sample(.N)], aes(PC1, PC2, colour = dataset)) +   # shuffled, so no dataset is drawn on top
    geom_point(size = 0.35, alpha = 0.6, stroke = 0) +
    scale_colour_manual(values = ds_pal, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = sprintf("PC1 (%.0f%%)", 100 * d$var_PC1[1]), y = sprintf("PC2 (%.0f%%)", 100 * d$var_PC2[1]),
         subtitle = c(uncorrected = "Uncorrected", corrected = "Batch-corrected")[cr]) +
    th + theme(axis.text = element_blank(), axis.ticks = element_blank())
}
row_a <- pca_panel("uncorrected") + labs(tag = "a") + pca_panel("corrected") + labs(tag = "b") +
  plot_layout(guides = "collect") & theme(legend.position = "right")

# c-d model types
mc  <- fread("results/tables/model_comparison_cv.csv")
mcs <- fread("results/tables/model_comparison_summary.csv")
m_lab <- c(glmnet = "Elastic net", ranger = "Random forest", knn = "k-nearest neighbours",
           svmRadial = "SVM (radial kernel)", pls = "PLS regression", xgboost = "Gradient boosting")
lev <- rev(m_lab[mcs$method])
mc  <- merge(mc, mcs[, .(method, kept)], by = "method")
mc[, method := factor(m_lab[method], levels = lev)]
mcs[, method := factor(m_lab[method], levels = lev)]
kept_pal <- c(`TRUE` = "#2a78d6", `FALSE` = "#a3a29d")

p_mc <- ggplot(mc, aes(r2, method)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = grey) +
  geom_point(aes(colour = kept), size = 0.9, alpha = 0.6, stroke = 0,
             position = position_jitter(width = 0, height = 0.15, seed = 1)) +
  geom_errorbar(data = mcs, aes(xmin = mean_r2, xmax = mean_r2, y = method), inherit.aes = FALSE,
                orientation = "y", width = 0.55, linewidth = 0.6, colour = ink) +
  scale_colour_manual(values = kept_pal, guide = "none") +
  labs(x = expression("Cross-validated" ~ italic(R)^2), y = NULL) + th
p_best <- ggplot(mcs, aes(n_best, method, fill = kept)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = kept_pal, labels = c(`TRUE` = "Kept", `FALSE` = "Dropped"), name = NULL,
                    breaks = c("TRUE", "FALSE")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 6, 2)) +
  labs(x = "Times ranked best (of 17)", y = NULL) +
  th + theme(axis.text.y = element_blank(), legend.position = "right")
row_b <- p_mc + labs(tag = "c") + p_best + labs(tag = "d") + plot_layout(widths = c(2.6, 1))

# e-f CV R-squared of the best model per drug and dataset (batch-corrected; identical uncorrected,
# since the correction only shifts each dataset by a constant)
grid <- CJ(feature_set = names(fs_lab), drug = drugs, dataset = datasets)
hm <- merge(grid, cv[best == TRUE & correction == "corrected", .(feature_set, drug, dataset, r2, method)],
            by = c("feature_set", "drug", "dataset"), all.x = TRUE)
short <- c(glmnet = "EN", svmRadial = "SVM", pls = "PLS")
hm[, `:=`(drug = factor(drug_lab[drug], levels = rev(drug_lab)), dataset = factor(ds_lab[dataset], levels = ds_lab[datasets]))]
lim <- ceiling(max(abs(hm$r2), na.rm = TRUE) * 10) / 10

hm_panel <- function(fs) {
  ggplot(hm[feature_set == fs], aes(dataset, drug, fill = r2)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(is.na(r2), "n.t.", sprintf("%.2f", r2))), size = 2.2, family = "Arial",
              colour = ink, nudge_y = 0.13) +
    geom_text(aes(label = ifelse(is.na(r2), "", short[method])), size = 1.8, family = "Arial",
              colour = "#5d5c58", nudge_y = -0.2) +
    scale_fill_gradient2(low = "#e34948", mid = "#f0efec", high = "#2a78d6", midpoint = 0,
                         limits = c(-lim, lim), na.value = "#f7f7f5", name = expression(italic(R)^2)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, subtitle = fs_lab[fs]) +
    th + theme(axis.line = element_blank(), axis.ticks = element_blank(),
               legend.key.height = unit(7, "mm"), legend.key.width = unit(2.5, "mm"))
}
row_c <- hm_panel("multiomics") + labs(tag = "e") + hm_panel("panel") + labs(tag = "f") + plot_layout(guides = "collect")

# Rows laid out separately (wrap_elements), so the long model names in c don't push a-b to the right
p1 <- wrap_elements(full = row_a) / wrap_elements(full = row_b) / wrap_elements(full = row_c) +
  plot_layout(heights = c(0.95, 0.8, 1.05))
save_fig(p1, "fig1_data_and_cv", 180, 160)

# ---- Figure 2: transfer between cell line datasets ----
cross <- cross[correction == "corrected"]
cross[, `:=`(drug = factor(drug_lab[drug], levels = drug_lab), feature_set = factor(fs_lab[feature_set], levels = fs_lab))]
dodge <- position_dodge(width = 0.7)
p2a <- ggplot(cross, aes(drug, r2, colour = feature_set)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = grey) +
  geom_point(size = 0.9, alpha = 0.55, stroke = 0, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7, seed = 1)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "errorbar", width = 0.22, linewidth = 0.7,
               position = dodge, show.legend = FALSE) +
  scale_colour_manual(values = fs_pal, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2, alpha = 1))) +
  labs(x = NULL, y = expression(italic(R)^2 ~ "on another dataset")) +
  th + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0))

# Mean R-squared over drugs for each train -> test pair, per feature set
pair <- cross[, .(r2 = mean(r2), n_drugs = .N), by = .(feature_set, train, test)]
pair <- merge(CJ(feature_set = unique(pair$feature_set), train = datasets, test = datasets), pair, all.x = TRUE)
pair[, `:=`(train = factor(ds_lab[train], levels = rev(ds_lab[datasets])), test = factor(ds_lab[test], levels = ds_lab[datasets]))]
pair_panel <- function(fs) {
  ggplot(pair[feature_set == fs], aes(test, train, fill = r2)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(is.na(r2), "", sprintf("%.2f", r2))), size = 2.4, family = "Arial", colour = ink) +
    scale_fill_gradient(low = "#f0efec", high = "#2a78d6", limits = c(0, 0.45), oob = scales::squish, na.value = "#ffffff",
                        name = expression("Mean" ~ italic(R)^2)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Test dataset", y = "Training dataset", subtitle = fs) +
    th + theme(axis.line = element_blank(), axis.ticks = element_blank(),
               axis.text.x = element_text(angle = 30, hjust = 1),
               legend.key.height = unit(6, "mm"), legend.key.width = unit(2.5, "mm"))
}
p2 <- p2a / (pair_panel("Multi-omics") + pair_panel("Gene panel") + plot_spacer() +
               plot_layout(guides = "collect", widths = c(1, 1, 0.35))) +
  plot_layout(heights = c(1.15, 1)) + tags
save_fig(p2, "fig2_cross_dataset", 180, 125)

# ---- Figure 3: PDX ----
ens <- pdx[model == "ensemble" & correction == "corrected"]
ens[, feature_set := factor(fs_lab[feature_set], levels = fs_lab)]
n_lab <- unique(ens[, .(drug, n)])[, setNames(sprintf("%s\n(n = %d)", drug_lab[drug], n), drug)]
ens[, drug := factor(n_lab[drug], levels = n_lab[drugs])]
bar_panel <- function(metric, ref, ylab, lims) {
  ggplot(ens, aes(drug, get(metric), fill = feature_set)) +
    geom_col(position = position_dodge(width = 0.62), width = 0.55, colour = "white", linewidth = 0.2) +
    geom_hline(yintercept = ref, linewidth = 0.3, colour = ink) +
    scale_fill_manual(values = fs_pal, name = NULL) +
    coord_cartesian(ylim = lims) +
    labs(x = NULL, y = ylab) + th + theme(axis.text.x = element_text(size = 6))
}
p3a <- bar_panel("spearman", 0, expression("Spearman" ~ rho), c(-0.3, 0.6))
p3b <- bar_panel("roc_auc_pd", 0.5, "ROC AUC (PD vs not PD)", c(0.3, 0.95)) +
  scale_y_continuous(breaks = seq(0.3, 0.9, 0.1))

# Predicted cell line AUC of the multi-omics ensemble, by mRECIST class
pp <- fread("results/tables/pdx_predictions.csv", na.strings = c("", "NA"))[
  feature_set == "multiomics" & correction == "corrected" & !is.na(mRECIST)]
pp[, class := factor(ifelse(mRECIST == "PD", "PD", "Not PD"), levels = c("Not PD", "PD"))]
auc <- pdx[model == "ensemble" & feature_set == "multiomics" & correction == "corrected", .(drug, roc_auc_pd)]
pp[, drug := factor(drug_lab[drug], levels = drug_lab)]
auc[, drug := factor(drug_lab[drug], levels = drug_lab)]
p3c <- ggplot(pp, aes(class, ensemble)) +
  geom_boxplot(aes(fill = class), width = 0.55, outlier.shape = NA, linewidth = 0.3, colour = ink, alpha = 0.9) +
  geom_point(size = 0.5, alpha = 0.55, stroke = 0, colour = ink, position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = auc, aes(x = 1.5, y = Inf, label = sprintf("ROC AUC %.2f", roc_auc_pd)), vjust = 1.3,
            size = 1.9, family = "Arial", colour = ink, inherit.aes = FALSE) +
  facet_wrap(~drug, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("Not PD" = "#f0efec", PD = "#86b6ef"), guide = "none") +
  labs(x = NULL, y = "Predicted AUC") + th +
  theme(axis.text = element_text(size = 5.5), axis.title = element_text(size = 7), strip.text = element_text(size = 7))

# d-f cancer type (06b): which types each drug was tested on, how much of the predictions and of
# the response it explains, Spearman over all PDX vs within cancer types. Multi-omics ensemble.
type_lab <- c(BRCA = "Breast", CM = "Melanoma", CRC = "Colorectal", NSCLC = "Lung (NSCLC)", PDAC = "Pancreatic",
              GC = "Gastric", unknown = "Unknown")
type_pal <- setNames(c("#d14f8f", "#1f1f1e", "#eb6834", "#2a78d6", "#1baf7a", "#e0b43a", "#c9c8c3"), type_lab)
ct_n <- fread("results/tables/pdx_cancer_type_counts.csv")
ct_s <- fread("results/tables/pdx_cancer_type_stats.csv")[feature_set == "multiomics" & correction == "corrected"]
ct_n[, `:=`(drug = factor(drug_lab[drug], levels = rev(drug_lab)), cancer_type = factor(type_lab[cancer_type], levels = type_lab))]
ct_s[, drug := factor(drug_lab[drug], levels = rev(drug_lab))]
pair_pal <- c("#bdbcb7", "#1c5cab")
dumbbell <- function(d, x, xlab) {
  ggplot(d, aes({{ x }}, drug)) +
    geom_line(aes(group = drug), linewidth = 0.5, colour = "#c9c8c3") +
    geom_point(aes(fill = variable), shape = 21, size = 2, stroke = 0.3, colour = "white") +
    scale_fill_manual(values = setNames(pair_pal, levels(d$variable)), name = NULL) +
    labs(x = xlab, y = NULL) + th + theme(legend.position = "bottom", legend.direction = "vertical")
}
p3d <- ggplot(ct_n, aes(n, drug, fill = cancer_type)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = type_pal, name = NULL, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = "PDX samples", y = NULL) + th + theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
ve <- melt(ct_s[n_types > 1], id.vars = "drug", measure.vars = c("var_obs", "var_pred"))
ve[, variable := factor(variable, labels = c("Observed PDX response", "Predicted AUC"))]
p3e <- dumbbell(ve, 100 * value, "Variance explained by cancer type (%)") +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 20))
vf <- melt(ct_s, id.vars = "drug", measure.vars = c("spearman_all", "spearman_within"))
vf[, variable := factor(variable, labels = c("All PDX", "Within cancer type"))]
p3f <- dumbbell(vf, value, expression("Spearman" ~ rho)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = grey)
p3f$layers <- c(p3f$layers[3], p3f$layers[1:2])   # zero line behind the points

# a and b side by side, one legend on top; c; d-f in one row
p3 <- (p3a + p3b + plot_layout(guides = "collect") & theme(legend.position = "top")) / p3c /
  (p3d + p3e + p3f + plot_layout(widths = c(1.2, 1, 1))) +
  plot_layout(heights = c(1, 1, 1.05)) + tags
save_fig(p3, "fig3_pdx", 180, 175)

# ---- Figure 4: top features ----
n_show <- 10
tf <- fread("results/tables/top_features.csv")[rank <= n_show]
omic_lab <- c(gex = "expr", cnv = "CNV", mut = "mut")
tf[, `:=`(drug = factor(drug_lab[drug], levels = drug_lab), key = sprintf("%s_%02d", drug, rank),
          direction = factor(ifelse(mean_weight > 0, "Resistance", "Sensitivity"), levels = c("Sensitivity", "Resistance")))]
# Gene symbols in italics (plotmath), Ensembl IDs without a symbol in plain text
tf[, label := ifelse(grepl("^ENSG", gene), sprintf('"%s"~"%s"', gene, omic_lab[omic]),
                     sprintf('italic("%s")~"%s"', gene, omic_lab[omic]))]
labs_map <- setNames(tf$label, tf$key)
tf[, key := factor(key, levels = rev(sort(unique(key))))]

p4 <- ggplot(tf, aes(1000 * mean_weight, key, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = ink) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d/%d", n_datasets, n_models), hjust = ifelse(mean_weight > 0, -0.25, 1.25)),
            size = 2, family = "Arial", colour = "#5d5c58") +
  facet_wrap(~drug, ncol = 3, scales = "free") +
  scale_y_discrete(labels = function(k) parse(text = labs_map[k])) +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  scale_fill_manual(values = c(Sensitivity = "#2a78d6", Resistance = "#e34948"), name = NULL,
                    labels = c(Sensitivity = "Higher value: more sensitive", Resistance = "Higher value: more resistant")) +
  labs(x = expression("Standardized weight (" * Delta * "predicted AUC per SD, " %*% 10^-3 * ")"), y = NULL) +
  th + theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5), panel.spacing.y = unit(4, "mm"))
save_fig(p4, "fig4_top_features", 180, 140)

print(summary)

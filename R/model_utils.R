# Helpers shared by the modeling scripts.

# Response table for one drug in one dataset (sample_id, value)
get_response <- function(data, drug_name) {
  resp <- data$drug[data$drug$drug == drug_name, ]
  resp[order(resp$sample_id), c("sample_id", "value")]
}

# Keep only some samples in every omic
subset_samples <- function(data, ids) {
  list(gex = data$gex[ids, , drop = FALSE],
       cnv = data$cnv[ids, , drop = FALSE],
       mut = data$mut[ids, , drop = FALSE])
}

# Regression metrics. R-squared here is 1 - SSE/SST (how close the values are),
# not caret's default squared correlation (only how well they line up).
# It can be negative when predictions are worse than just guessing the mean.
reg_metrics <- function(obs, pred) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]; pred <- pred[ok]
  same <- sd(pred) == 0   # constant predictions have no correlation
  data.frame(
    n        = length(obs),
    r2       = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    pearson  = if (same) NA_real_ else cor(obs, pred),
    spearman = if (same) NA_real_ else cor(obs, pred, method = "spearman"),
    rmse     = sqrt(mean((obs - pred)^2))
  )
}

# Train one model with caret. No preProcess: glmnet standardizes the features
# internally (using the training data only) and random forest doesn't need it.
train_model <- function(x, y, method, ctrl, n_threads = 1) {
  if (method == "glmnet") {
    caret::train(x, y, method = "glmnet", trControl = ctrl, tuneLength = 5)
  } else if (method == "ranger") {
    mtry <- unique(round(sqrt(ncol(x)) * c(0.5, 1, 2)))
    caret::train(x, y, method = "ranger", trControl = ctrl,
                 num.trees = 500, num.threads = n_threads,
                 tuneGrid = expand.grid(mtry = mtry, splitrule = "variance",
                                        min.node.size = 5))
  } else {
    stop("unknown method: ", method)
  }
}

# Predict from a saved model on any dataset (features filtered the same way as in training)
predict_saved <- function(model, data, ids) {
  x <- apply_feature_filter(subset_samples(data, ids), model$filter)
  as.numeric(predict(model$fit, newdata = x))
}

model_path <- function(drug, dataset, method) {
  file.path("models", "task1", paste0(drug, "_", dataset, "_", method, ".rds"))
}

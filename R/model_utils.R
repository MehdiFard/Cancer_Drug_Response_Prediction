# Helpers shared by the modeling scripts.

# Response table for one drug in one dataset (sample_id, value)
get_response <- function(data, drug_name) {
  resp <- data$drug[data$drug$drug == drug_name, ]
  resp[order(resp$sample_id), c("sample_id", "value")]
}

# Keep only some samples in every omic (the panel data has no gex)
subset_samples <- function(data, ids) {
  omics <- intersect(c("gex", "cnv", "mut"), names(data))
  lapply(setNames(omics, omics), function(o) data[[o]][ids, , drop = FALSE])
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

# ROC AUC via Mann-Whitney: probability that a positive sample scores higher than a negative one.
# 0.5 = no better than chance, 1 = perfect separation.
roc_auc <- function(score, positive) {
  ok <- !is.na(score) & !is.na(positive)
  score <- score[ok]; positive <- positive[ok]
  n1 <- sum(positive); n0 <- sum(!positive)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score)[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# caret's built-in xgboost models break with xgboost 3.x, so this is a small
# custom caret model that calls xgboost directly. The booster is wrapped in a list
# because caret tries to add fields to it, which xgboost 3 doesn't allow.
xgb_caret <- list(
  label = "xgboost", library = "xgboost", type = "Regression",
  parameters = data.frame(parameter = c("nrounds", "max_depth", "eta"), class = "numeric",
                          label = c("Rounds", "Max depth", "Learning rate")),
  grid = function(x, y, len = NULL, search = "grid") {
    expand.grid(nrounds = 150, max_depth = c(2, 4), eta = 0.1)
  },
  fit = function(x, y, wts, param, lev, last, classProbs, nthread = 1, ...) {
    params <- list(objective = "reg:squarederror", max_depth = param$max_depth, eta = param$eta,
                   subsample = 0.8, colsample_bytree = 0.2, tree_method = "hist", nthread = nthread)
    booster <- xgboost::xgb.train(params = params, nrounds = param$nrounds, verbose = 0,
                                  data = xgboost::xgb.DMatrix(as.matrix(x), label = y))
    list(booster = booster)
  },
  predict = function(modelFit, newdata, submodels = NULL) {
    predict(modelFit$booster, xgboost::xgb.DMatrix(as.matrix(newdata)))
  },
  prob = NULL,
  sort = function(x) x[order(x$max_depth), ]
)

# Train one model with caret.
# Scaling: glmnet standardizes the features itself and the tree models don't need it.
# kNN, SVM and PLS get center + scale inside train, so it's learned on each training
# fold only ("zv" first drops columns that happen to be constant in a fold).
train_model <- function(x, y, method, ctrl, n_threads = 1) {
  scale_pp <- c("zv", "center", "scale")
  switch(method,
    glmnet = caret::train(x, y, method = "glmnet", trControl = ctrl, tuneLength = 5),
    ranger = caret::train(x, y, method = "ranger", trControl = ctrl,
                          num.trees = 500, num.threads = n_threads,
                          tuneGrid = expand.grid(mtry = unique(round(sqrt(ncol(x)) * c(0.5, 1, 2))),
                                                 splitrule = "variance", min.node.size = 5)),
    knn = caret::train(x, y, method = "knn", trControl = ctrl, preProcess = scale_pp,
                       tuneGrid = data.frame(k = c(5, 10, 20, 40))),
    svmRadial = caret::train(x, y, method = "svmRadial", trControl = ctrl, preProcess = scale_pp,
                             tuneLength = 5),
    pls = caret::train(x, y, method = "simpls", trControl = ctrl, preProcess = scale_pp,
                       tuneGrid = data.frame(ncomp = 1:10)),
    xgboost = caret::train(x, y, method = xgb_caret, trControl = ctrl, nthread = n_threads),
    stop("unknown method: ", method)
  )
}

# Methods that use their own threads; the others run CV folds in parallel instead
threaded_methods <- c("ranger", "xgboost")

# Drop what a trained model doesn't need for predicting, before saving it.
# - the training matrix (caret keeps it, glmnet keeps another copy in its call)
# - PLS keeps the training data and a formula table of features x features (>1 GB here)
#   and needs them to predict, so we predict straight from its coefficients instead:
#   y = (x - mean of x) %*% B + mean of y
slim_fit <- function(fit) {
  fit$trainingData <- NULL
  if (is.list(fit$finalModel)) fit$finalModel$call <- NULL   # SVM is an S4 object without one
  if (inherits(fit$finalModel, "mvr")) {
    fit$finalModel[c("model", "terms", "fitted.values", "residuals")] <- NULL
    fit$modelInfo$predict <- function(modelFit, newdata, submodels = NULL) {
      B <- modelFit$coefficients[, 1, modelFit$ncomp]
      as.vector(sweep(as.matrix(newdata), 2, modelFit$Xmeans) %*% B + modelFit$Ymeans)
    }
  }
  fit
}

# Predict from a saved model on any dataset (features filtered the same way as in training)
predict_saved <- function(model, data, ids) {
  x <- apply_feature_filter(subset_samples(data, ids), model$filter)
  as.numeric(predict(model$fit, newdata = x))
}

# One folder per feature set and correction, e.g. models/multiomics_uncorrected/
model_path <- function(drug, dataset, method, feature_set, correction) {
  file.path("models", paste0(feature_set, "_", correction), paste0(drug, "_", dataset, "_", method, ".rds"))
}

# Processed data of one dataset for one feature set ("multiomics" or "panel"), batch
# "uncorrected" (from 01, <dataset>_clean_<feature set>.rds) or "corrected" (from 02)
load_data <- function(dataset, feature_set, correction = "uncorrected") {
  stage <- c(uncorrected = "clean", corrected = "corrected")[[correction]]
  readRDS(file.path("data/processed", paste0(dataset, "_", stage, "_", feature_set, ".rds")))
}

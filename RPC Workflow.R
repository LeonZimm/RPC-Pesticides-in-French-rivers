# RPC - Addtioninal Pesticide Predicitons in French Rivers using large Monoroting Data
setwd("C:/Users/leonz/Documents/Ecotox_Studium/4. Semester/RPC/RPC Pesticides in French Rivers")

# required packages
library(tidyverse)
library(caret)
library(data.table)
library(ranger)
library(foreach)
library(doParallel)
library(MLmetrics)
library(precrec)
library(mltools)

# Load Data
sa10b400_train <- readRDS("Data/sa10b400_train.rds")
sa10b400_test <- readRDS("Data/sa10b400_test.rds")
sa10b400_validation <- readRDS("Data/sa10b400_validation.rds")

# Hyperparameter Optimization Process
hp_tun_fct <- function(data,
                       param_name,
                       param_values,
                       fixed_params = list(),
                       kfolds = 5,
                       seed = 25) {
  
  # Table for results
  res_log <- data.table()
  # Create folds
  set.seed(seed)
  folds <- createFolds(data$det, k = kfolds, list = TRUE, returnTrain = FALSE)
  
  for (val in param_values) {
    cat("Tuning", param_name, "=", val,
        "| Time:", format(Sys.time(), "%H:%M:%S"),
        "\n")
    
    
    for (k in 1:kfolds) {
      start <- proc.time()[[3]]
      cat(" Fold", k,
          " | Time:", format(Sys.time(), "%H:%M:%S"),
          "\n")
      
      # Prepare data
      idx <- setdiff(seq_len(nrow(data)), folds[[k]])
      train <- data[idx]
      test <- data[folds[[k]]]
      
      # Build dynamic parameter list
      param_list <- list(
        formula = det ~.,
        data = train,
        respect.unordered.factors = "order",
        probability = TRUE,
        num.threads = detectCores() - 1,
        seed = seed
      )
      
      # Add tuning parameter
      param_list[[param_name]] <- val
      
      # Add fixed parameters
      param_list <- c(param_list, fixed_params)
      
      # Train model
      model <- do.call(ranger, param_list)
      
      # Predict, calibrate probabilities (isotonic regression), compute AUC/PR AUC and threshold 
      res_list <- evaluation_fct(model, test)
      
      end <- proc.time()[[3]]
      # Error log
      res_log <- rbind(res_log, data.table(param_value = val,
                                           auc = res_list$auc,
                                           pr_auc = res_list$pr_auc,
                                           f1 = res_list$f1,
                                           mcc =  res_list$mcc,
                                           time = (end - start)/60))
    }
    
    
  }
  
  setnames(res_log, "param_value", param_name)
  
  res <- res_log[, .(auc = mean(auc),
                     pr_auc = mean(pr_auc),
                     f1 = mean(f1),
                     mcc = mean(mcc),
                     time = mean(time)), 
                 by = param_name]
  
  
  if (class(param_values) == "numeric") {
    if (param_name == "min.node.size") {
      tunval_plot <- ndsz_num_hypa_fct(res, param_name)
    } else{
      tunval_plot <- num_hypa_fct(res, param_name)  
    }
  } else {
    tunval_plot <- char_hypa_fct(res, param_name)
  }
  
  output <- list("res_log" = res_log[], "res" = res, "tunval" = tunval_plot[[1]], "plot" = tunval_plot[[2]])
  
  return(output)
}


# Hyperparameter Value Selection
num_hypa_fct <- function(res, para){
  best <- max(res$mcc)
  best_zone <- max(best) - 0.01 * max(best)
  
  if (!any((res$mcc > (best_zone)) & (res[[para]] < res[[para]][which.max(res$mcc)]))) {
    para_tun = res[[para]][which.max(res$mcc)]
  } else {
    para_tun = res[[para]][res[[para]] == min(res[[para]][(res$mcc > (best_zone)) & (res[[para]] < res[[para]][which.max(res$mcc)])])]
  }
  
  plot <- ggplot(res, aes(x = !!sym(para), y = mcc, group = 1)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = best_zone, color = "#cc0000")
  
  output <- list("para_tun" = para_tun, "plot" = plot)
  return(output)
}
ndsz_num_hypa_fct <- function(res, para){
  best <- max(res$mcc)
  best_zone <- max(best) - 0.01 * max(best)
  
  if (!any((res$mcc > (best_zone)) & (res[[para]] > res[[para]][which.max(res$mcc)]))) {
    para_tun = res[[para]][which.max(res$mcc)]
  } else {
    para_tun = res[[para]][res[[para]] == max(res[[para]][(res$mcc > (best_zone)) & (res[[para]] > res[[para]][which.max(res$mcc)])])]
  }
  
  plot <- ggplot(res, aes(x = !!sym(para), y = mcc, group = 1)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = best_zone, color = "#cc0000")
  
  output <- list("para_tun" = para_tun, "plot" = plot)
  return(output)
}
char_hypa_fct <- function(res, para) {
  best <- max(res$mcc)
  best_zone <- max(best) - 0.01 * max(best)
  
  if (sum(res$mcc > best_zone) == 1) {
    para_tun = res[[para]][which.max(res$mcc)]
  } else {
    para_tun = res[[para]][which.min(res$mcc[res$mcc > best_zone])]
  }
  
  plot <- ggplot(res, aes(x = !!sym(para), y = mcc, group = 1)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = best_zone, color = "#cc0000")
  
  output <- list("para_tun" = para_tun, "plot" = plot)
  return(output)
}

# Predictions and Performance Evaluation
evaluation_fct <- function(model, test) {
  
  # Prediction
  prob <- predict(model, test)$predictions
  
  # AUC
  mm <- mmdata(prob[, 1],
               test$det,
               posclass = "D")
  evalmod <- evalmod(mm)
  auc <- auc(evalmod)[1, 4]
  pr_auc <- auc(evalmod)[2, 4]
  # F1
  pred <- as.factor(ifelse(prob[, 1] >= 0.5, "D", "ND"))
  f1 <- F1_Score(y_pred = pred,
                 y_true = test$det,
                 positive = "D")
  # MCC
  mcc <- mltools::mcc(as.factor(pred),
                      test$det)
  
  list("auc" = auc,
       "pr_auc" = pr_auc,
       "f1" = f1, 
       "mcc" = mcc,
       "pred" = pred,
       "evalmod" = as.data.table(evalmod))
}
IR_eval_fct <- function(model, validation, test, n_cores) {
  
  # Prediction
  prob_val <- predict(model, validation)$predictions[, "D"]
  
  thresholds <- seq(0.01, 0.99, 0.02)
  
  # Set up parallelization
  n_cores <- n_cores
  clus <- makeCluster(n_cores)
  on.exit(stopCluster(clus), add = TRUE)
  registerDoParallel(clus)
  
  thre_res <- foreach(t = thresholds,
                      .combine = rbind,
                      .packages = c("MLmetrics", "mltools", "data.table")) %dopar% {
                        
                        pred <- ifelse(prob_val >= t, "D", "ND")
                        
                        # MCC
                        mcc <- mltools::mcc(as.factor(pred),
                                            validation$det)
                        
                        data.table(threshold = t, mcc = mcc)
                      }
  
  best_t <- thre_res[which.max(mcc), threshold]
  
  prob_test <- predict(model, test)$predictions[, "D"]
  pred <- ifelse(prob_test >= best_t, "D", "ND")
  pred_num <- ifelse(pred == "D", 1, 0)
  test_num <- ifelse(test$det == "D", 1, 0)
  # AUC
  mm <- mmdata(pred_num,
               test_num,
               posclass = 1)
  evalmod <- evalmod(mm)
  
  auc <- auc(evalmod)[1, 4]
  pr_auc <- auc(evalmod)[2, 4]
  # F1
  f1 <- F1_Score(y_pred = pred,
                 y_true = test$det,
                 positive = "D")
  # MCC
  mcc <- mltools::mcc(as.factor(pred),
                      test$det)
  
  list("auc" = auc,
       "pr_auc" = pr_auc, 
       "f1" = f1,
       "mcc" = mcc,
       "threshold" = best_t,
       "pred" = pred,
       "evalmod" = as.data.table(evalmod))
}

# 1. Feature Selection using Recursive Feature Elimination -------------------------------------------------------

rfs_data <- copy(sa10b400_train)            

# Variables
n_trees <- 120                     
kfolds <- 5                    
t_max <- ncol(rfs_data) - 1        


# Store results for analysis
rfs_eval_log <- data.table()
rfs_feat_log <- data.table()

for (t in 1:t_max) {
  cat("Variables left", ncol(rfs_data), 
      "| Iteration", t,
      "| Time:", format(Sys.time(), "%H:%M:%S"),
      "\n")
  
  # Store k-fold variable importance measurements
  impu_list <- list()
  perm_list <- list()
  
  # Create k-folds 
  set.seed(25)
  folds_rfs <- createFolds(rfs_data$det, 
                           k = kfolds, 
                           list = TRUE, 
                           returnTrain = FALSE)
  
  start <- proc.time()[[3]]
  for (k in 1:kfolds) {
    cat(" Fold", k,
        " | Time:", format(Sys.time(), "%H:%M:%S"),
        "\n")
    
    # Prepare data
    rfs_index <- folds_rfs[[k]]
    rfs_train <- rfs_data[-rfs_index]
    rfs_test <- rfs_data[rfs_index]
    
    # RF for impurity importance
    rfs_impu <- ranger(det ~ .,
                       data = rfs_train,
                       num.trees = n_trees,
                       importance = "impurity",
                       respect.unordered.factors = "order",
                       num.threads = detectCores() - 1,
                       seed = 25)
    
    # RF for permutation importance
    rfs_perm <- ranger(det ~ .,
                       data = rfs_train,
                       num.trees = n_trees,
                       importance = "permutation",
                       respect.unordered.factors = "order",
                       num.threads = detectCores() - 1,
                       probability = TRUE,
                       seed = 25)
    
    # Predict, calibrate probabilities (isotonic regression), compute AUC/PR AUC and threshold
    res_list <- evaluation_fct(rfs_perm, rfs_test)
    
    # Variable importance table + error log
    impu_list[[k]] <- as.data.table(rfs_impu$variable.importance, keep.rownames = TRUE)
    perm_list[[k]] <- as.data.table(rfs_perm$variable.importance, keep.rownames = TRUE)
    
    rfs_eval_log <- rbind(rfs_eval_log, data.table(iter = t,
                                                   auc = res_list$auc,
                                                   pr_auc = res_list$pr_auc,
                                                   f1 = res_list$f1,
                                                   mcc = res_list$mcc))
  }
  end <- proc.time()[[3]]
  time_rfs <- end - start                           
  
  # Calculate mean of cross validated importance measurements
  impu_dt <- rbindlist(impu_list)
  perm_dt <- rbindlist(perm_list)
  
  impu_mean <- impu_dt[, .(impu = mean(V2)), by = V1]
  perm_mean <- perm_dt[, .(perm = mean(V2)), by = V1]
  
  va_impo_dt <- merge(impu_mean, perm_mean, by = "V1") 
  
  # Normalizing importance measurements and scoring features
  z_score <- function(x) (x - mean(x))/sd(x)
  
  va_impo_dt[, impu_z := z_score(impu)]
  va_impo_dt[, perm_z := z_score(perm)]
  va_impo_dt[, score := 0.5*impu_z + 0.5*perm_z] 
  va_impo_dt <- va_impo_dt[order(score)]
  
  # Removing 5% of the worst features
  feat_rem <- va_impo_dt$V1[1:(0.05*nrow(va_impo_dt))]
  rfs_data <- rfs_data[, (feat_rem) := NULL]
  
  # Storing results
  rfs_feat_log <- rbind(rfs_feat_log, data.table(iter = t,
                                                 n_features = nrow(va_impo_dt),
                                                 feature_rm = feat_rem,
                                                 time = time_rfs))
  
  if (ncol(rfs_data) <= 1) {
    break
  }
  
  print(t)
}

rfs_res <- rfs_eval_log[, .(auc = mean(auc),
                            pr_auc = mean(pr_auc),
                            f1 = mean(f1),
                            mcc = mean(mcc)),
                        by = iter]

rfs_res <- merge(rfs_res, rfs_feat_log, by = "iter", all.x = TRUE)

rfs_plot <- ggplot(rfs_res, aes(as.factor(n_features), mcc)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = max(rfs_res$mcc) - 0.01 * max(rfs_res$mcc), color = "#cc0000")

# Identify the columns to remove
cols_to_remove <- rfs_res$feature_rm[1:(which.max(rfs_res$mcc) - 1)]
# Get remaining column names
rfs_features <- setdiff(colnames(sa10b400_train), cols_to_remove)

rfs <- list("eval_log" = rfs_eval_log,
            "feat_log" = rfs_feat_log,
            "res" = rfs_res,
            "plot" = rfs_plot,
            "remaining_features" = rfs_features)
saveRDS(rfs, "RPC_Full_Run/rfs.rds")


# 2. Hyperparameter Optimization using a Sequential Method ------------------------------------------------------

hp_tun_data <- copy(sa10b400_train[, ..rfs_features])


# 2.1 Number of trees -----------------------------------------------------

ntree_val <- c(20, 50, 100, 125, 200, 500)

ntree <- hp_tun_fct(hp_tun_data,
                    "num.trees",
                    ntree_val,
                    fixed_params = list())
saveRDS(ntree, "RPC_Full_Run/ntree.rds")

# 2.2 Mtry ----------------------------------------------------------------

mtry_val <- seq(2, ncol(hp_tun_data) - 1, by = 2)

mtry <- hp_tun_fct(hp_tun_data,
                   "mtry",
                   mtry_val,
                   fixed_params = list(num.trees = ntree$tunval))
saveRDS(mtry, "RPC_Full_Run/mtry.rds")


# 2.3 Sample fraction ---------------------------------------------------------

safr_val <- seq(1, 0.5, by = -0.1)

safr <- hp_tun_fct(hp_tun_data,
                   "sample.fraction",
                   safr_val,
                   fixed_params = list(num.trees = ntree$tunval,
                                       mtry = mtry$tunval))
saveRDS(safr, "RPC_Full_Run/safr.rds")


# 2.4 Minimal node size ---------------------------------------------------

ndsz_val <- seq(1, 6, 1)

ndsz <- hp_tun_fct(hp_tun_data,
                   "min.node.size",
                   ndsz_val,
                   fixed_params = list(num.trees = ntree$tunval,
                                       mtry = mtry$tunval,
                                       sample.fraction = safr$tunval,
                                       replace = repl$tunval))
saveRDS(ndsz, "RPC_Full_Run/ndsz.rds")


# 2.5 Splitting rule ------------------------------------------------------

sprl_val <- c("gini", "extratrees", "hellinger")

sprl <- hp_tun_fct(hp_tun_data,
                   "splitrule",
                   sprl_val,
                   fixed_params = list(num.trees = ntree$tunval,
                                       mtry = mtry$tunval,
                                       sample.fraction = safr$tunval,
                                       replace = repl$tunval,
                                       min.node.size = ndsz$tunval))
saveRDS(sprl, "RPC_Full_Run/sprl.rds")


# 3. Null-Model --------------------------------------------------

basli_model_ranger <- ranger(det ~ .,
                             data = sa10b400_train,
                             num.trees = 120,
                             respect.unordered.factors = "order",
                             importance = "permutation",
                             num.threads = detectCores() - 1,
                             probability = TRUE,
                             seed = 25)

basli_res_list <- evaluation_fct(basli_model_ranger, sa10b400_validation)
cf_basli <- confusionMatrix(sa10b400_validation$det,
                            as.factor(basli_res_list$pred),
                            positive = "D")

basli_model <- list("res" = basli_res_list,
                    "confu_mat" = cf_basli,
                    "model" = basli_model_ranger)
saveRDS(basli_model, "RPC_Results/basli_model.rds")


# 4. Optimized model --------------------------------------------------

opt_data <- copy(sa10b400_train[, ..rfs_features])

opt_model_ranger <- ranger(det ~ .,
                           data = opt_data,
                           num.trees = 50,
                           mtry = 8,
                           sample.fraction = 0.9,
                           min.node.size = 6,
                           splitrule = "extratrees",
                           importance = "permutation",
                           probability = TRUE,
                           respect.unordered.factors = "order",
                           num.threads = detectCores() - 1,
                           seed = 25)

opt_res_list <- IR_eval_fct(opt_model_ranger, sa10b400_test, sa10b400_validation, 5)
cf_opt <- confusionMatrix(sa10b400_validation$det,
                          as.factor(opt_res_list$pred),
                          positive = "D")

opt_model <- list("res" = opt_res_list,
                  "confu_mat" = cf_opt,
                  "model" = opt_model_ranger,
                  "features" = rfs_features)
saveRDS(opt_model, "RPC_Full_Run/opt_model.rds")

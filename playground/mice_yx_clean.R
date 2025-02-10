## this is for imputation method: mice_yy ##
# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))

# config for this imputation method
model <- "CS"
model_type <- "LM"
horizon <- "7 days"
imputation <- "mice_yx"
# path_data_miss <- data_path_play_dyn_miss
# selected variables
vars_selected <- vars_impute_before

preds_name <- paste("preds", model, horizon, imputation, sep = "_")
results_name <- paste("results", model, horizon, imputation, sep = "_")
coefs_name <- paste("coefs", model, horizon, imputation, sep = "_")
timings_name <- paste("timings", model, horizon, imputation, sep = "_")

# get filenames for datasets with mising values
datasets_files <- list.files(data_path_play_dyn_miss, 
                             recursive = TRUE, full.names = TRUE)
train_files <- datasets_files[str_detect(datasets_files, "train")]
test_files <- datasets_files[str_detect(datasets_files, "test")]

# keep results
predictions <- init_preds_BASE_DYN()
results <- init_results_BASE_DYN()
coefs <- init_coefs()
timings <- init_timings()

# set.seed(2024)

# impute train and test
for (f in train_files){
  
  print(f)
  start_time <- Sys.time()
  
  # impute train
  load(f) # loads train data named data_train
  test_file <- str_replace(f, "train", "test") # corresponding test set file
  load(test_file) 
  
  
  # vars_select + vars_id_outcome
  data_train <-  data_train %>% select(vars_selected, vars_id_outcome, c("LAB_WBC_count_last","LAB_WBC_Neutrophils_last"))
  data_test <-  data_test %>% select(vars_selected, vars_id_outcome, c("LAB_WBC_count_last","LAB_WBC_Neutrophils_last"))
  
  # binarize
  cat_cols <- "MS_medical_specialty"
  
  bin_model <- data_train %>% 
    make_col_binary_drop(cat_cols, dropped_levels = list(MS_medical_specialty = "Other"))
  data_train <- bin_model$data
  bin_model <- data_test %>% 
    make_col_binary_drop(cat_cols, dropped_levels = list(MS_medical_specialty = "Other"))
  data_test <- bin_model$data
  
  
  
  # train data set
  
  data_train$LAB_urea_last_log <- log(data_train$LAB_urea_last)  
  data_train$LAB_creatinine_last_log <- log(data_train$LAB_creatinine_last)
  data_train$LAB_bilirubin_last_log <- log(data_train$LAB_bilirubin_last)
  data_train$LAB_ciclosporin_last_log <- log(data_train$LAB_ciclosporin_last)
  
  # test data set
  
  data_test$LAB_urea_last_log <- log(data_test$LAB_urea_last)
  data_test$LAB_creatinine_last_log <- log(data_test$LAB_creatinine_last)
  data_test$LAB_bilirubin_last_log <- log(data_test$LAB_bilirubin_last)
  data_test$LAB_ciclosporin_last_log <- log(data_test$LAB_ciclosporin_last)
  
  
  #first setup predictorMatrix so that Y is not used to impute
  data_train <- data_train %>% 
    filter(LM <= 30) %>% 
    select(vars_mice3, 
           c("functioneelDossierNr", "CAT_catheter_episode", "LM", "type", "eventtime"), 
           c("LAB_WBC_count_last","LAB_WBC_Neutrophils_last")) %>%
    select(-c("LAB_is_neutropenia"))
  
  data_train$MS_medical_specialty_bin_drop_Cardiac <- as.factor(data_train$MS_medical_specialty_bin_drop_Cardiac)
  data_train$MS_medical_specialty_bin_drop_Traumatology <- as.factor(data_train$MS_medical_specialty_bin_drop_Traumatology)
  data_train$MS_medical_specialty_bin_drop_Pediatrics <- as.factor(data_train$MS_medical_specialty_bin_drop_Pediatrics)
  data_train$ADM_admission_source_binary_all_Home <- as.factor(data_train$ADM_admission_source_binary_all_Home)
  
  set.seed(2024)
  start_time_imputation <- Sys.time()
  dummyMice <- mice(data_train, m = 1, maxit = 0)
  predMat <- dummyMice$predictorMatrix
  predMat[c("functioneelDossierNr"), ] <- 0
  predMat[,c("functioneelDossierNr")] <- 0
  
  imputeMethod <- make.method(data_train)
  logreg_vars <- c("ADM_admission_source_binary_all_Home", "MS_medical_specialty_bin_drop_Cardiac", "MS_medical_specialty_bin_drop_Traumatology", "MS_medical_specialty_bin_drop_Pediatrics")
  imputeMethod[logreg_vars] <- "logreg"  # Logistic regression for binary variables

  
  MI <- mice(data_train, m = 10, predictorMatrix = predMat, method = imputeMethod, maxit = 5, print = F) 
  # mi_x_train <- mice::complete(imp.train, action = "long", include = TRUE)
  # imputed.train <- as.mids(mi_x_train)
  data_test <- data_test %>% 
    filter(LM <= 30) %>% 
    select(vars_mice3, 
           c("functioneelDossierNr", "CAT_catheter_episode", "LM", "type", "eventtime"), 
           c("LAB_WBC_count_last","LAB_WBC_Neutrophils_last")) %>%
    select(-c("LAB_is_neutropenia"))
  
  data_test$MS_medical_specialty_bin_drop_Cardiac <- as.factor(data_test$MS_medical_specialty_bin_drop_Cardiac)
  data_test$MS_medical_specialty_bin_drop_Traumatology <- as.factor(data_test$MS_medical_specialty_bin_drop_Traumatology)
  data_test$MS_medical_specialty_bin_drop_Pediatrics <- as.factor(data_test$MS_medical_specialty_bin_drop_Pediatrics)
  data_test$ADM_admission_source_binary_all_Home <- as.factor(data_test$ADM_admission_source_binary_all_Home)
  
  data_test_x <- data_test %>% 
    mutate(type = ifelse(!is.na(type), NA, type),
           eventtime = ifelse(!is.na(eventtime), NA, eventtime))
  
  data_test_y <- data_test %>% 
    # filter(LM <= 30) %>% 
    select(c("functioneelDossierNr", "CAT_catheter_episode", "LM", "type", "eventtime")) %>%
    filter(eventtime > LM) 
  
  imp.test <- mice.mids(MI, newdata = data_test_x, maxit = 5)
  
 
  imp_train_list <- list()
  
  for (i in 1:10) {
    
    complete_data <- create_new_mice_yx2(complete(MI, action = i), FALSE)  # Extract complete dataset
    
    imp_train_list[[i]] <- complete_data
    
  }
  
  
  imp_test_list <- list()
  
  for (i in 1:10) {
    
    complete_data <- create_new_mice_yx2(complete(imp.test, action = i), data_test_y)  # Extract complete dataset
    
    imp_test_list[[i]] <- complete_data
    
  }
  time_imputation <- as.numeric(difftime(Sys.time(), start_time_imputation, units = "secs"))
  
  # specify the formula
  form <- paste0(" ~ cluster(id) + ",paste0(predictors_col, collapse = " + "))
  
  # List to store the models
  model_list <- list()
  
  model_info_list <- list()
  
  # Fit model on each imputed dataset
  start_time_run_model <- Sys.time()
  for (i in 1:length(imp_train_list)) {
    mod <- riskRegression::CSC(update.formula(prodlim::Hist(eventtime,type,Tstart)~., form), imp_train_list[[i]], method = "breslow")
    model_list[[i]] <- mod
    model_info_list[[i]] <- list(
      "coefficients" = stats::coef(mod),
      "baseline_hazards" = lapply(mod$models, function(model) survival::basehaz(model, centered = FALSE)),
      "model_terms" = lapply(mod$models, function(model) model[["terms"]])
    )
  }
  time_run_model <- as.numeric(difftime(Sys.time(), start_time_run_model, units = "secs"))
  
  preds_list <- list()
  
  start_time_predict <- Sys.time()
  for (i in 1:length(imp_train_list)) {
    preds <- predictLM_CSC(model_info_list[[i]], predictors_col, imp_test_list[[i]], horizon = 7, primary_event = 1)
    preds_list[[i]] <- preds
  }
  time_predict <- as.numeric(difftime(Sys.time(), start_time_predict, units = "secs"))
  
  
  data_train_bin <- data_train %>%
    mutate(type = if_else((type == "CLABSI"|type == "Death"|type == "Discharge") & eventtime > LM + 7, "Censored", type),
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
           Tstart = LM,
           eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)) %>%
    filter(!eventtime <= LM) %>% 
    mutate(type = if_else(type == "Censored", 0, if_else(type == "CLABSI", 1, if_else(type == "Death", 2, 3))))
  
  
  
  data_test_bin <- data_test %>%
    mutate(type = if_else((type == "CLABSI"|type == "Death"|type == "Discharge") & eventtime > LM + 7, "Censored", type),
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
           Tstart = LM,
           eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)) %>%
    filter(!eventtime <= LM) %>% 
    mutate(type = if_else(type == "Censored", 0, if_else(type == "CLABSI", 1, if_else(type == "Death", 2, 3))))
  
  
  data_test_bin$pred_test <- rowMeans(do.call(cbind, preds_list))
  
  coefs_list <- list()
  for (i in 1:length(model_info_list)) {
    coef <- model_info_list[[i]]$coefficients$`Cause 1`
    coefs_list[[i]] <- coef
  }
  
  coef_cs_model <- data.frame(rowMeans(do.call(cbind, coefs_list)))
  colnames(coef_cs_model) <-"value"
  
  # observed risk within 7 days y_true_cat (categorical: "CLABSI", "no_CLABSI", "Discharge", "Death", "Censored")
  # observed risk within 7 days y_test (binary: 0/1)
  data_test_bin$y_true_cat <- if_else(data_test_bin$type == 1, "CLABSI", if_else(data_test_bin$type == 2, "Death", if_else(data_test_bin$type == 3, "Discharge", "Censored")))
  
  data_test_bin$y_test <- if_else(data_test_bin$type == 1, 1, 0)
  
  
  predictions <- predictions %>% 
    add_row(preds = data_test_bin$pred_test,
            y_true_cat = data_test_bin$y_true_cat,
            y_true_time = data_test_bin$eventtime,
            train_set = str_replace(f, data_path_play_dyn_miss, ""),
            test_set = str_replace(test_file, data_path_play_dyn_miss, ""),
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = data_test_bin$LM,
            functioneelDossierNr = data_test_bin$functioneelDossierNr,
            CAT_catheter_episode = data_test_bin$CAT_catheter_episode)
  
  # performance metrics of the model
  
  metrics <- data_test_bin %>% group_by(LM) %>% summarise(
    AUROC = c(c_statistic(y_test, pred_test)),
    slope = calibration_slope(y_test, pred_test),
    intercept = calibration_large(y_test, pred_test),
    OE_ratio = oe_ratio(y_test, pred_test),
    ECE = ECE(y_test, pred_test),
    ECI = ECI(y_test, pred_test),
    BS = brier_score(y_test, pred_test),
    Scaled_BS = scaled_brier_score(y_test, pred_test)) %>%
    pivot_longer(cols= AUROC:Scaled_BS,
                 names_to="metric",
                 values_to="value")
  
  results <- results %>%
    add_row(train_set = str_replace(f, data_path_play_dyn_miss, ""),
            test_set = str_replace(test_file, data_path_play_dyn_miss, ""),
            metric = metrics$metric,
            value = metrics$value,
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = metrics$LM)
  
  coefs <- coefs %>% 
    add_row(variable = rownames(coef_cs_model),
            value = coef_cs_model$value,
            train_set = str_replace(f, data_path_play_dyn_miss, ""),
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = NA)
  
  timings <- timings %>% 
    add_row(type = c("impute", "build model", "predict"),
            value = c(time_imputation, time_run_model, time_predict),
            train_set = str_replace(f, data_path_play_dyn_miss, ""),
            model = imputation,
            LM = NA_real_) 
  
  message(sprintf("DONE in %s minutes.", 
                  difftime(Sys.time(), start_time, units = "mins") %>% as.numeric()))
}

# save predictions 
save(predictions, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/predictions/", preds_name))
save(results, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/performances/", results_name))
save(coefs, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/coefs/", coefs_name))
save(timings, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/timings/", timings_name))

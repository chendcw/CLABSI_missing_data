## this is for imputation method: single mice ##
# source all files in the R directory
files_to_source <- list.files("/data/leuven/342/vsc34229/clabsi_r/R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))

# config for this imputation method
model <- "CS"
model_type <- "LM"
horizon <- "7 days"
imputation <- "single"
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

# train_files <- train_files[1:10]
# test_files <- test_files[1:10]

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
  
  
  data_train <- data_train %>% 
    filter(LM <= 30) %>% 
    select(vars_mice3, 
           c("functioneelDossierNr", "CAT_catheter_episode", "LM", "type", "eventtime"), 
           c("LAB_WBC_count_last","LAB_WBC_Neutrophils_last")) %>%
    select(-c("LAB_is_neutropenia"))
  
  # data_train_bin$MS_medical_specialty_bin_drop_Cardiac <- as.factor(data_train_bin$MS_medical_specialty_bin_drop_Cardiac)
  # data_train_bin$MS_medical_specialty_bin_drop_Traumatology <- as.factor(data_train_bin$MS_medical_specialty_bin_drop_Traumatology)
  # data_train_bin$MS_medical_specialty_bin_drop_Pediatrics <- as.factor(data_train_bin$MS_medical_specialty_bin_drop_Pediatrics)
  # data_train$LAB_is_neutropenia <- as.factor(data_train$LAB_is_neutropenia)
  data_train$ADM_admission_source_binary_all_Home <- as.factor(data_train$ADM_admission_source_binary_all_Home)
  # data_train$MS_medical_specialty <- as.factor(data_train$MS_medical_specialty)
  data_train$MS_medical_specialty_bin_drop_Cardiac <- as.factor(data_train$MS_medical_specialty_bin_drop_Cardiac)
  data_train$MS_medical_specialty_bin_drop_Traumatology <- as.factor(data_train$MS_medical_specialty_bin_drop_Traumatology)
  data_train$MS_medical_specialty_bin_drop_Pediatrics <- as.factor(data_train$MS_medical_specialty_bin_drop_Pediatrics)
  
  
  # data_train$CAT_lumens_CVC <- as.factor(data_train$CAT_lumens_CVC)
  # data_train$CAT_lumens_Tunneled_CVC <- as.factor(data_train$CAT_lumens_Tunneled_CVC)
  # data_train$CAT_lumens_PICC <- as.factor(data_train$CAT_lumens_PICC)
  # data_train$CARE_NEU_GCS_score_last <- as.factor(data_train$CARE_NEU_GCS_score_last)
  
  data_test <- data_test %>% 
    filter(LM <= 30) %>% 
    select(vars_mice3, 
           c("functioneelDossierNr", "CAT_catheter_episode", "LM", "type", "eventtime"), 
           c("LAB_WBC_count_last","LAB_WBC_Neutrophils_last")) %>%
    select(-c("LAB_is_neutropenia"))
  
  # data_test_bin$MS_medical_specialty_bin_drop_Cardiac <- as.factor(data_test_bin$MS_medical_specialty_bin_drop_Cardiac)
  # data_test_bin$MS_medical_specialty_bin_drop_Traumatology <- as.factor(data_test_bin$MS_medical_specialty_bin_drop_Traumatology)
  # data_test_bin$MS_medical_specialty_bin_drop_Pediatrics <- as.factor(data_test_bin$MS_medical_specialty_bin_drop_Pediatrics)
  # data_test$LAB_is_neutropenia <- as.factor(data_test$LAB_is_neutropenia)
  data_test$ADM_admission_source_binary_all_Home <- as.factor(data_test$ADM_admission_source_binary_all_Home)
  # data_test$MS_medical_specialty <- as.factor(data_test$MS_medical_specialty)
  data_test$MS_medical_specialty_bin_drop_Cardiac <- as.factor(data_test$MS_medical_specialty_bin_drop_Cardiac)
  data_test$MS_medical_specialty_bin_drop_Traumatology <- as.factor(data_test$MS_medical_specialty_bin_drop_Traumatology)
  data_test$MS_medical_specialty_bin_drop_Pediatrics <- as.factor(data_test$MS_medical_specialty_bin_drop_Pediatrics)
  
  # data_test$CAT_lumens_CVC <- as.factor(data_test$CAT_lumens_CVC)
  # data_test$CAT_lumens_Tunneled_CVC <- as.factor(data_test$CAT_lumens_Tunneled_CVC)
  # data_test$CAT_lumens_PICC <- as.factor(data_test$CAT_lumens_PICC)
  # data_test$CARE_NEU_GCS_score_last <- as.factor(data_test$CARE_NEU_GCS_score_last)
  
  set.seed(2024)
  start_time_imputation <- Sys.time()
  dummyMice <- mice(data_train, m = 1, maxit = 0)
  predMat <- dummyMice$predictorMatrix
  predMat[c("eventtime", "type", "functioneelDossierNr"), ] <- 0
  predMat[,c("eventtime", "type", "functioneelDossierNr")] <- 0
  imputeMethod <- make.method(data_train)
  # polyreg_vars <- c("CAT_lumens_CVC", "CAT_lumens_Tunneled_CVC", "CAT_lumens_PICC", "CARE_NEU_GCS_score_last", "MS_medical_specialty")
  # logreg_vars <- c("ADM_admission_source_binary_all_Home")
  logreg_vars <- c("ADM_admission_source_binary_all_Home", "MS_medical_specialty_bin_drop_Cardiac", "MS_medical_specialty_bin_drop_Traumatology", "MS_medical_specialty_bin_drop_Pediatrics")
  
  # imputeMethod[polyreg_vars] <- "polyreg"     # polyreg regression for count variables
  imputeMethod[logreg_vars] <- "logreg"  # Logistic regression for binary variables
  # imputeMethod[linear_vars] <- "pmm" 
  ## do not allow Y in imputation model 
  MI <- mice(data_train, m = 1, predictorMatrix = predMat, method = imputeMethod, maxit = 1, print = F) 
  
  
  imp.test <- mice.mids(MI, newdata = data_test, maxit = 1)
  
  data_train_imputed <- complete(MI)  # Extract complete dataset
  
  data_train_imputed$LM_1 <- data_train_imputed$LM/30
  data_train_imputed$LM_2 <- (data_train_imputed$LM/30)^2
  data_train_imputed$MS_is_ICU_unit_1 <-  data_train_imputed$MS_is_ICU_unit * data_train_imputed$LM_1
  data_train_imputed$MS_is_ICU_unit_2 <-  data_train_imputed$MS_is_ICU_unit * data_train_imputed$LM_2
  
  data_train_imputed <- data_train_imputed %>%
    mutate(CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath = ifelse(CAT_catheter_type_binary_all_Tunneled_CVC==1|CAT_catheter_type_binary_all_Port_a_cath==1, 1, 0))
  
  data_train_imputed <- data_train_imputed %>%
    mutate(CAT_lumens_Total = CAT_lumens_CVC + 
             CAT_lumens_Tunneled_CVC +
             CAT_lumens_Dialysis_CVC +
             CAT_lumens_Port_a_cath +
             CAT_lumens_PICC) 
  
  # data_train_imputed <- data_train_imputed %>%
  #   mutate(CAT_lumens_Total = as.numeric(as.character(CAT_lumens_CVC)) + 
  #            as.numeric(as.character(CAT_lumens_Tunneled_CVC)) +
  #            as.numeric(as.character(CAT_lumens_Dialysis_CVC)) +
  #            as.numeric(as.character(CAT_lumens_Port_a_cath)) +
  #            as.numeric(as.character(CAT_lumens_PICC))) 
  
  data_train_imputed <- data_train_imputed %>% 
    mutate(LAB_is_neutropenia = as.numeric(LAB_WBC_count_last * 1000 < 500 |
                                             LAB_WBC_Neutrophils_last * 1000 < 500)) 
  
  
  data_test_imputed <- complete(imp.test)
  
  data_test_imputed$LM_1 <- data_test_imputed$LM/30
  data_test_imputed$LM_2 <- (data_test_imputed$LM/30)^2
  data_test_imputed$MS_is_ICU_unit_1 <-  data_test_imputed$MS_is_ICU_unit * data_test_imputed$LM_1
  data_test_imputed$MS_is_ICU_unit_2 <-  data_test_imputed$MS_is_ICU_unit * data_test_imputed$LM_2
  
  data_test_imputed <- data_test_imputed %>%
    mutate(CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath = ifelse(CAT_catheter_type_binary_all_Tunneled_CVC==1|CAT_catheter_type_binary_all_Port_a_cath==1, 1, 0))
  
  data_test_imputed <- data_test_imputed %>%
    mutate(CAT_lumens_Total = CAT_lumens_CVC + 
             CAT_lumens_Tunneled_CVC +
             CAT_lumens_Dialysis_CVC +
             CAT_lumens_Port_a_cath +
             CAT_lumens_PICC) 
  
  # data_test_imputed <- data_test_imputed %>%
  #   mutate(CAT_lumens_Total = as.numeric(as.character(CAT_lumens_CVC)) + 
  #            as.numeric(as.character(CAT_lumens_Tunneled_CVC)) +
  #            as.numeric(as.character(CAT_lumens_Dialysis_CVC)) +
  #            as.numeric(as.character(CAT_lumens_Port_a_cath)) +
  #            as.numeric(as.character(CAT_lumens_PICC))) 
  
  data_test_imputed <- data_test_imputed %>% 
    mutate(LAB_is_neutropenia = as.numeric(LAB_WBC_count_last * 1000 < 500 |
                                             LAB_WBC_Neutrophils_last * 1000 < 500)) 
  time_imputation <- as.numeric(difftime(Sys.time(), start_time_imputation, units = "secs"))
  
  # binarize
  # cat_cols <- "MS_medical_specialty"
  
  # bin_model <- data_train_imputed %>% 
  #   make_col_binary_drop(cat_cols, dropped_levels = list(MS_medical_specialty = "Other"))
  # data_train_bin <- bin_model$data
  # data_test_bin <- data_test_imputed %>% 
  #   make_col_binary_drop(cat_cols, dropped_levels = list(MS_medical_specialty = "Other"))
  # data_test_bin <- data_test_bin$data
  
  # outcome: Surv(eventtime, type)
  
  data_train_bin <- data_train_imputed %>%
    mutate(type = if_else((type == "CLABSI"|type == "Death"|type == "Discharge") & eventtime > LM + 7, "Censored", type),
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
           Tstart = LM,
           eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)) %>%
    filter(!eventtime <= LM) %>% 
    mutate(type = if_else(type == "Censored", 0, if_else(type == "CLABSI", 1, if_else(type == "Death", 2, 3))))
  
  
  
  data_test_bin <- data_test_imputed %>%
    mutate(type = if_else((type == "CLABSI"|type == "Death"|type == "Discharge") & eventtime > LM + 7, "Censored", type),
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
           Tstart = LM,
           eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)) %>%
    filter(!eventtime <= LM) %>% 
    mutate(type = if_else(type == "Censored", 0, if_else(type == "CLABSI", 1, if_else(type == "Death", 2, 3))))
  
  # data_train_bin$LAB_is_neutropenia <- as.numeric(as.character(data_train_bin$LAB_is_neutropenia))
  data_train_bin$ADM_admission_source_binary_all_Home <- as.numeric(as.character(data_train_bin$ADM_admission_source_binary_all_Home))
  # data_train_bin$CARE_NEU_GCS_score_last <- as.numeric(as.character(data_train_bin$CARE_NEU_GCS_score_last))
  
  
  # data_test_bin$LAB_is_neutropenia <- as.numeric(as.character(data_test_bin$LAB_is_neutropenia))
  data_test_bin$ADM_admission_source_binary_all_Home <- as.numeric(as.character(data_test_bin$ADM_admission_source_binary_all_Home))
  # data_test_bin$CARE_NEU_GCS_score_last <- as.numeric(as.character(data_test_bin$CARE_NEU_GCS_score_last))
  
  # specify the formula
  form <- paste0(" ~ cluster(id) + ",paste0(predictors_col, collapse = " + "))
  
  
  # fit landmark supermodel
  # simple
  # use LMsupercs0
  start_time_run_model <- Sys.time()
  LMsupercs0 <- riskRegression::CSC(update.formula(prodlim::Hist(eventtime,type,Tstart)~., form), data_train_bin, method = "breslow")
  time_run_model <- as.numeric(difftime(Sys.time(), start_time_run_model, units = "secs"))
  
  coef_cs_model <- coef(LMsupercs0$models$`Cause 1`)
  
  # predict on test set
  start_time_predict <- Sys.time()
  model_info <- list(
    "coefficients" = stats::coef(LMsupercs0),
    "baseline_hazards" = lapply(LMsupercs0$models, function(mod) survival::basehaz(mod, centered = FALSE)),
    "model_terms" = lapply(LMsupercs0$models, function(mod) mod[["terms"]])
  )
  
  # rows since 80708 have problem
  data_test_bin$pred_test <- predictLM_CSC(model_info, predictors_col, newdata = data_test_bin, horizon = 7, primary_event = 1)
  time_predict <- as.numeric(difftime(Sys.time(), start_time_predict, units = "secs"))
  
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
    add_row(variable = names(coef_cs_model),
            value = coef_cs_model,
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
save(predictions, file = paste0("/data/leuven/342/vsc34229/clabsi_r/playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/predictions/", preds_name))
save(results, file = paste0("/data/leuven/342/vsc34229/clabsi_r/playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/performances/", results_name))
save(coefs, file = paste0("/data/leuven/342/vsc34229/clabsi_r/playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/coefs/", coefs_name))
save(timings, file = paste0("/data/leuven/342/vsc34229/clabsi_r/playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/timings/", timings_name))

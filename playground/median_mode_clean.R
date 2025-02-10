## this is for imputation method: median/mode ##
# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))

# config for this imputation method
model <- "CS"
model_type <- "LM"
horizon <- "7 days"
imputation <- "median_mode"
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


# impute train and test
for (f in train_files){
  
  print(f)
  start_time <- Sys.time()
  
  # impute train
  load(f) # loads train data named data_train
  test_file <- str_replace(f, "train", "test") # corresponding test set file
  load(test_file) 
  
  # vars_select + vars_id_outcome
  data_train <-  data_train %>% select(vars_selected, vars_id_outcome)
  data_test <-  data_test %>% select(vars_selected, vars_id_outcome)
  # Impute missing values using median/mode per LM
  
  # classify columns into vars_impute_mode & vars_impute_median
  
  # int [0/1] --> mode
  # num [0/1] --> mode
  # int [not 0/1] --> mode
  # num [not 0/1] --> median
  # cat/chr --> mode
  start_time_imputation <- Sys.time()
  vars <- classify_variables(data_train)
  vars_impute_median <- vars$vars_impute_median
  vars_impute_median <- setdiff(vars_impute_median, "LM")
  vars_impute_mode <- vars$vars_impute_mode
  vars_impute_mode <- setdiff(vars_impute_mode, "LM")
  
  # CAT_lumens_Dialysis_CVC & CAT_lumens_Port_a_cath are exceptional
  vars_impute_median <- setdiff(vars_impute_median,c("CAT_lumens_Dialysis_CVC","functioneelDossierNr","CAT_catheter_episode","eventtime"))
  vars_impute_mode <- setdiff(vars_impute_mode, c("CAT_lumens_Port_a_cath","type"))
  
  
  # obtain the median/mode value in training set
  # save the median/mode values into a dataframe
  median_values <- data_train %>%
    filter(LM <= 30) %>%
    group_by(LM) %>%
    summarise(across(all_of(vars_impute_median), median, na.rm = TRUE)) 
  
  mode_values <- data_train %>%
    filter(LM <= 30) %>%
    group_by(LM) %>%
    summarise(across(all_of(vars_impute_mode), Mode)) 
  
  median_mode_values <- median_values %>%
    left_join(mode_values, by = "LM")
  
  # take only the data until LM 30
  
  data_train_imputed <- data_train %>%
    filter(LM <= 30) %>%
    group_by(LM) %>%
    mutate_at(vars(vars_impute_median), funs(impute_median)) %>%
    mutate_at(vars(vars_impute_mode), funs(impute_mode)) %>%
    ungroup()
  
  # apply the saved median/mode values to the test set for imputation
  
  data_test_imputed <- data_test %>% left_join(median_mode_values, by = "LM") %>% filter(LM <= 30)
  
  
  # replace NA values in test dataset with corresponding median values
  for (var in vars_impute_median) {
    data_test_imputed[[paste0(var, ".x")]] <- ifelse(is.na(data_test_imputed[[paste0(var, ".x")]]), data_test_imputed[[paste0(var, ".y")]], data_test_imputed[[paste0(var, ".x")]])
  }
  
  
  # replace NA values in test dataset with corresponding mode values
  for (var in vars_impute_mode) {
    data_test_imputed[[paste0(var, ".x")]] <- ifelse(is.na(data_test_imputed[[paste0(var, ".x")]]), data_test_imputed[[paste0(var, ".y")]], data_test_imputed[[paste0(var, ".x")]])
  }
  time_imputation <- as.numeric(difftime(Sys.time(), start_time_imputation, units = "secs"))
  
  
  # remove unnecessary columns
  data_test_imputed <- select(data_test_imputed, -ends_with(".y"))
  
  
  data_test_imputed <- data_test_imputed %>%
    rename_with(~gsub("\\.x$", "", .), ends_with(".x"))

  # binarize
  cat_cols <- "MS_medical_specialty"
  
  bin_model <- data_train_imputed %>% 
    make_col_binary_drop(cat_cols, dropped_levels = list(MS_medical_specialty = "Other"))
  data_train_bin <- bin_model$data
  data_test_bin <- data_test_imputed %>% 
    make_col_binary_drop(cat_cols, dropped_levels = list(MS_medical_specialty = "Other"))
  data_test_bin <- data_test_bin$data
  
  
  # outcome: Surv(eventtime, type)
  
  data_train_bin <- data_train_bin %>%
    mutate(type = if_else((type == "CLABSI"|type == "Death"|type == "Discharge") & eventtime > LM + 7, "Censored", type),
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
           Tstart = LM,
           eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)) %>%
    filter(!eventtime <= LM) %>% 
    mutate(type = if_else(type == "Censored", 0, if_else(type == "CLABSI", 1, if_else(type == "Death", 2, 3))))
  
  
  
  data_test_bin <- data_test_bin %>%
    mutate(type = if_else((type == "CLABSI"|type == "Death"|type == "Discharge") & eventtime > LM + 7, "Censored", type),
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
           Tstart = LM,
           eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)) %>%
    filter(!eventtime <= LM) %>% 
    mutate(type = if_else(type == "Censored", 0, if_else(type == "CLABSI", 1, if_else(type == "Death", 2, 3))))
  
  # train data set
  
  data_train_bin <- create_new_variable(data_train_bin) %>% select(predictors_col, vars_id_outcome, c("id", "Tstart"))
  
  # test data set
  
  data_test_bin <- create_new_variable(data_test_bin) %>% select(predictors_col, vars_id_outcome, c("id", "Tstart"))
 

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
save(predictions, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/predictions/", preds_name))
save(results, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/performances/", results_name))
save(coefs, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/coefs/", coefs_name))
save(timings, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/timings/", timings_name))



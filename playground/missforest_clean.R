## this is for imputation method: missforest ##
# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))

# config for this imputation method
model <- "CS"
model_type <- "LM"
horizon <- "7 days"
imputation <- "missforest"
# path_data_miss <- data_path_play_dyn_miss
# selected variables
vars_selected <- vars_impute_before

preds_name <- paste("preds", model, horizon, imputation, sep = "_")
results_name <- paste("results", model, horizon, imputation, sep = "_")
coefs_name <- paste("coefs", model, horizon, imputation, sep = "_")
timings_name <- paste("timings", model, horizon, imputation, sep = "_")
OOB_errs_name <- paste("OOB_errs", model, horizon, imputation, sep = "_")

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
OOB_errs <- tibble(iteration = integer(),
                   variable = character(),
                   MSE = double(),
                   NMSE = double(),
                   MER = double(),
                   macro_F1 = double(),
                   F1_score = double(),
                   train_set = character())


# impute train and test
for (f in train_files){
  
  print(f)
  start_time <- Sys.time()
  
  # impute train
  load(f) # loads train data named data_train
  test_file <- str_replace(f, "train", "test") # corresponding test set file
  load(test_file) 
  
  data_train <- data_train %>% 
    select(all_of(c(vars_impute_selected, vars_id_outcome))) %>% 
    arrange(functioneelDossierNr, CAT_catheter_episode, LM)
  
  data_test <- data_test %>% 
    select(all_of(c(vars_impute_selected, vars_id_outcome))) %>% 
    arrange(functioneelDossierNr, CAT_catheter_episode, LM)
  
  
  # keep outcome separately
  data_train_outcome <- data_train %>% 
    select(eventtime, type, functioneelDossierNr, CAT_catheter_episode, LM)
  
  data_test_outcome <- data_test %>% 
    select(eventtime, type, functioneelDossierNr, CAT_catheter_episode, LM)
  
  data_train <- data_train %>% 
    select(-c(eventtime, type))
  
  data_train <- data_train %>% 
    arrange(functioneelDossierNr, CAT_catheter_episode, LM)
  
  data_test <- data_test %>% 
    select(-c(eventtime, type))
  
  data_test <- data_test %>% 
    arrange(functioneelDossierNr, CAT_catheter_episode, LM)
  
  (table_miss <- map(data_train, ~sum(is.na(.))) %>% 
      as_tibble() %>% 
      pivot_longer(cols = everything(), names_to = "feature", values_to = "number_missing") %>% 
      mutate(number_present = nrow(data_train) - number_missing,
             percentage_missing = number_missing/dim(data_train)[[1]],
             percentage_missing = round(percentage_missing, 4)) %>% 
      arrange(desc(number_missing)) %>% 
      as.data.frame())
  
  
  support_vars <- table_miss %>% 
    filter(percentage_missing == 0) %>% 
    pull(feature)
  support_vars <- support_vars[!support_vars %in% c("functioneelDossierNr", 
                                                    "type", "eventtime")]
  
  start_time_imputation <- Sys.time()
  predictor_matrix <- missForestPredict::create_predictor_matrix(data_train)
  predictor_matrix[support_vars,] <- 0
  
  # make integer columns (e.g: binary) double - bagging will impute with values between 0 and 1
  column_types <- sapply(data_train, class)
  integer_column_types <- names(column_types[column_types == "integer"])
  integer_columns <- integer_column_types[!integer_column_types %in% c("functioneelDossierNr", "CAT_catheter_episode", "LM")]
  
  data_train <- data_train %>% 
    mutate_at(all_of(integer_columns), as.double)
  
  data_test <- data_test %>% 
    mutate_at(all_of(integer_columns), as.double)
  
  # change verbose = TRUE to save OOB error and runtime
  missForest_model <- missForestPredict::missForest(data_train, verbose = FALSE,
                                           initialization = "median/mode",
                                           predictor_matrix = predictor_matrix,
                                           num.trees = 500,
                                           maxiter = 100)

  
  # save OOB err
  OOB_errs <- OOB_errs %>% 
    add_row(missForest_model$OOB_err %>% add_column(train_set = f))

  # missForest_model <- fit_imputation_model_missForest(data_train, verbose = FALSE, 
  #                                                     # proportion_usable_cases = c(0.95, 0.05),
  #                                                     predictor_matrix = predictor_matrix)
  
  # impute train
  data_train_imputed <- missForest_model$ximp
  
  # impute test data 
  data_test_imputed <- missForestPredict::missForestPredict(missForest_model, 
                                                            newdata = data_test)
  time_imputation <- as.numeric(difftime(Sys.time(), start_time_imputation, units = "secs"))
  
  # data_test_imputed <- predict_imputation_missForest(missForest_model, 
  #                                                    newdata = data_test, 
  #                                                    verbose = FALSE)
  
  # binarize
  cat_cols <- "MS_medical_specialty"
  
  bin_model <- data_train_imputed %>% 
    make_col_binary_drop(cat_cols, dropped_levels = list(MS_medical_specialty = "Other"))
  data_train_bin <- bin_model$data
  data_test_bin <- data_test_imputed %>% 
    make_col_binary_drop(cat_cols, dropped_levels = list(MS_medical_specialty = "Other"))
  data_test_bin <- data_test_bin$data
  
  # outcome: Surv(eventtime, type)
  
  data_train_bin <- cbind(data_train_bin, data_train_outcome %>% select(-c("LM", "CAT_catheter_episode", "functioneelDossierNr")))
  # data_train_bin <- cbind(data_train_bin, data_train_outcome %>% select(-c("LM", "CAT_catheter_episode")))
  
  data_test_bin <- cbind(data_test_bin, data_test_outcome %>% select(-c("LM", "CAT_catheter_episode", "functioneelDossierNr")))
  
  
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
  
  
  data_train_bin <- create_new_variable(data_train_bin) %>% filter(LM<=30) %>% select(predictors_col, vars_id_outcome, c("id", "Tstart"))
  data_test_bin <- create_new_variable(data_test_bin) %>% filter(LM<=30) %>% select(predictors_col, vars_id_outcome, c("id", "Tstart"))
  
  
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
save(OOB_errs, file = paste0("playground/3_Imputation_for_Dynamic_CLABSI_Prediction_in_EHR_Data/OOB_errs/", OOB_errs_name))


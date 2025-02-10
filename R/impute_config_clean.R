# imputation study configuration
library(mice)

data_path_play_dyn_miss <- "data_for_models/2012_2013/missing/DYN/"

# vars with identification & outcome
vars_id_outcome <- c("functioneelDossierNr", 
                     "CAT_catheter_episode",
                     # generate "id" using "functioneelDossierNr" & "CAT_catheter_episode", 
                     "LM", 
                     # generate "Tstart" using "LM"
                     "type", 
                     "eventtime")


# vars before logrithm & binarization & combination

vars_impute_before <- c("PAT_age",
                        "ADM_admission_source_binary_all_Home",
                        "CAT_catheter_type_binary_all_CVC",
                        "CAT_catheter_type_binary_all_Tunneled_CVC",
                        "CAT_catheter_type_binary_all_Port_a_cath",
                        "CAT_catheter_type_binary_all_PICC",
                        "CAT_lumens_CVC",
                        "CAT_lumens_Tunneled_CVC",
                        "CAT_lumens_Dialysis_CVC",
                        "CAT_lumens_Port_a_cath",
                        "CAT_lumens_PICC",
                        "MED_7d_TPN",
                        "MED_L2_7d_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE",
                        "MED_L2_7d_L01_ANTINEOPLASTIC_AGENTS",
                        "CARE_NEU_GCS_score_last",
                        "MS_medical_specialty",
                        "MS_is_ICU_unit",
                        "CARE_VS_temperature_max",
                        "LAB_is_neutropenia",
                        "LAB_urea_last",
                        "LAB_creatinine_last",
                        "LAB_bilirubin_last",
                        "LAB_ciclosporin_last",
                        "LAB_pH_last")



# vars after logrithm & binarization & combination


vars_impute_after <- c("PAT_age",
                       "ADM_admission_source_binary_all_Home",
                       "CAT_catheter_type_binary_all_CVC",
                       "CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath",
                       "CAT_catheter_type_binary_all_PICC",
                       "CAT_lumens_Total",
                       "MED_7d_TPN",
                       "MED_L2_7d_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE",
                       "MED_L2_7d_L01_ANTINEOPLASTIC_AGENTS",
                       "CARE_NEU_GCS_score_last",
                       "MS_medical_specialty_bin_drop_Cardiac",
                       "MS_medical_specialty_bin_drop_Traumatology",
                       "MS_medical_specialty_bin_drop_Pediatrics",
                       "MS_is_ICU_unit",
                       "CARE_VS_temperature_max",
                       "LAB_is_neutropenia",
                       "LAB_urea_last_log",
                       "LAB_creatinine_last_log",
                       "LAB_bilirubin_last_log",
                       "LAB_ciclosporin_last_log",
                       "LAB_pH_last")


# vars for building landmark cause-specific model
predictors_col <- c("PAT_age",
                    "ADM_admission_source_binary_all_Home",
                    "CAT_catheter_type_binary_all_CVC",                             
                    "CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath",                     
                    "CAT_catheter_type_binary_all_PICC", 
                    "CAT_lumens_Total",
                    "MED_7d_TPN",   
                    "MED_L2_7d_L01_ANTINEOPLASTIC_AGENTS",                          
                    "MED_L2_7d_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE", 
                    "CARE_NEU_GCS_score_last",
                    "MS_medical_specialty_bin_drop_Cardiac",
                    "MS_medical_specialty_bin_drop_Traumatology",
                    "MS_medical_specialty_bin_drop_Pediatrics",
                    "MS_is_ICU_unit",                                               
                    "CARE_VS_temperature_max",                                      
                    "LAB_is_neutropenia",
                    "LAB_urea_last_log",
                    "LAB_creatinine_last_log",
                    "LAB_bilirubin_last_log",
                    "LAB_ciclosporin_last_log",
                    "LAB_pH_last",
                    "LM_1",
                    "LM_2",
                    "MS_is_ICU_unit_1",
                    "MS_is_ICU_unit_2")




# for logrithm etc.
create_new_variable <- function(data) {
  # Here, you can define your logic to create the new variable
  # For demonstration purposes, let's create a new variable named 'new_var' with random values
  data$LAB_urea_last_log <- log(data$LAB_urea_last)  # Example: Generating random values
  data$LAB_creatinine_last_log <- log(data$LAB_creatinine_last)
  data$LAB_bilirubin_last_log <- log(data$LAB_bilirubin_last)
  data$LAB_ciclosporin_last_log <- log(data$LAB_ciclosporin_last)
  
  data$LM_1 <- data$LM/30
  data$LM_2 <- (data$LM/30)^2
  data$MS_is_ICU_unit_1 <-  data$MS_is_ICU_unit * data$LM_1
  data$MS_is_ICU_unit_2 <-  data$MS_is_ICU_unit * data$LM_2
  
  data <- data %>%
    mutate(CAT_lumens_Total = CAT_lumens_CVC + CAT_lumens_Tunneled_CVC + CAT_lumens_Dialysis_CVC + CAT_lumens_Port_a_cath + CAT_lumens_PICC) 
  
  data <- data %>%
    mutate(CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath = ifelse(CAT_catheter_type_binary_all_Tunneled_CVC==1|CAT_catheter_type_binary_all_Port_a_cath==1, 1, 0))
  
  return(data)
}



classify_variables <- function(data) {
  vars_impute_mode <- character(0)  # Initialize vector for variables to impute using mode
  vars_impute_median <- character(0)  # Initialize vector for variables to impute using median
  
  for (col in names(data)) {
    unique_vals <- unique(data[[col]])
    
    if (any(is.na(unique_vals))) {
      unique_vals <- unique_vals[!is.na(unique_vals)]
    }
    
    if (is.integer(data[[col]]) || is.numeric(data[[col]])) {
      if (all(unique_vals == 0 | unique_vals == 1)) {
        vars_impute_mode <- c(vars_impute_mode, col)
      } else {
        vars_impute_median <- c(vars_impute_median, col)
      }
    } else if (is.factor(data[[col]]) || is.character(data[[col]])) {
      vars_impute_mode <- c(vars_impute_mode, col)
    }
  }
  
  list(
    vars_impute_mode = vars_impute_mode,
    vars_impute_median = vars_impute_median
  )
}

impute_median <- function(x) {
  median_value <- median(x, na.rm = TRUE)
  ifelse(is.na(x), median_value, x)
}


# Function to impute missing values with mode for a single variable within each group
impute_mode <- function(x) {
  # Check the class of the input vector
  class_x <- class(x)
  
  if (all(is.na(x))) {
    # If all values are NA, return as is
    return(x)
  } else {
    # Compute the mode value
    mode_value <- names(sort(table(x), decreasing = TRUE))[1]
    
    # Impute missing values with mode value while preserving the class
    if (class_x == "integer") {
      return(ifelse(is.na(x), as.integer(mode_value), x))
    } else if (class_x == "numeric") {
      return(ifelse(is.na(x), as.numeric(mode_value), x))
    } else if (class_x == "character") {
      return(ifelse(is.na(x), as.character(mode_value), x))
    } else {
      return(x)  # Return as is for other classes
    }
  }
}


# Mode function
Mode <- function(x) {
  x <- x[!is.na(x)]  # Remove NA values
  ux <- unique(x)
  if (length(ux) == 0) return(NA)  # Return NA if all values are NA
  ux[which.max(tabulate(match(x, ux)))]
}

# for missing indicator approach, generate 1 for all values which are missing
fill_missing_values <- function(data) {
  for (i in 1:ncol(data)) {
    if (is.character(data[[i]])) {
      data[[i]][is.na(data[[i]])] <- "Empty"
    } else if (is.numeric(data[[i]])) {
      data[[i]][is.na(data[[i]])] <- as.numeric(99)
    } else if (is.integer(data[[i]])) {
      data[[i]][is.na(data[[i]])] <- as.integer(99)
    }
  }
  return(data)
}


# single stochastic regression

single_stochastic_regression_impute <- function(data, target_var, predictors) {
  # Filter rows with missing values for the target variable
  data_missing <- data %>% filter(is.na(.data[[target_var]]))
  data_nonmissing <- data %>% filter(!is.na(.data[[target_var]]))
  
  # Choose the imputation model based on the type of target variable
  if (is.numeric(data[[target_var]]) && 
      all(data[[target_var]] >= 0, na.rm = TRUE) &&
      all(data[[target_var]] == floor(data[[target_var]]), na.rm = TRUE) && # Ensure integers
      max(data[[target_var]], na.rm = TRUE) > 1) {
    
    # For count data, use Poisson regression
    # model <- glm.nb(as.formula(paste(target_var, "~", paste(predictors, collapse = " + "))),
    #                 data = data_nonmissing)
    model<- glm(as.formula(paste(target_var, "~", paste(predictors, collapse = " + "))),
                data = data_nonmissing, family = poisson(link = "log"))
    
    # Get predictions and add residual noise for stochasticity
    predictions <- predict(model, newdata = data_missing, type = "response")
    stochastic_imputations <- rpois(length(predictions), lambda = predictions)

    
  } else if (all(data[[target_var]] %in% c(0, 1, NA), na.rm = TRUE)) {
    # For binary data, use logistic regression
    model <- glm(as.formula(paste(target_var, "~", paste(predictors, collapse = " + "))),
                 data = data_nonmissing, family = binomial(link = "logit"))
    
    # Get predictions and add stochastic noise for stochasticity
    predictions <- predict(model, newdata = data_missing, type = "response")
    stochastic_imputations <- rbinom(length(predictions), size = 1, prob = predictions)
    
  } else {
    # For continuous data, use linear regression
    model <- lm(as.formula(paste(target_var, "~", paste(predictors, collapse = " + "))),
                data = data_nonmissing)
    
    # Get predictions and add residual noise for stochasticity
    residual_sd <- sd(residuals(model))
    predictions <- predict(model, newdata = data_missing)
    stochastic_imputations <- predictions + rnorm(length(predictions), mean = 0, sd = residual_sd)
  }
  
  # Replace missing values in the original data with stochastic imputations
  data[[target_var]][is.na(data[[target_var]])] <- stochastic_imputations
  
  return(data)
}


single_stochastic_regression_impute <- function(data_train, data_test, target_var, predictors) {
  # Split data into rows with missing and non-missing values for the target variable
  data_test_missing <- data_test %>% filter(is.na(.data[[target_var]]))
  data_train_nonmissing <- data_train %>% filter(!is.na(.data[[target_var]]))
  
  # Choose the imputation model based on the type of target variable
  if (is.numeric(data_train[[target_var]]) && 
      all(data_train[[target_var]] >= 0, na.rm = TRUE) &&
      all(data_train[[target_var]] == floor(data_train[[target_var]]), na.rm = TRUE) && # Ensure integers
      max(data_train[[target_var]], na.rm = TRUE) > 1) {
    
    # For count data, use Poisson regression
    model <- glm(as.formula(paste(target_var, "~", paste(predictors, collapse = " + "))),
                 data = data_train_nonmissing, family = poisson(link = "log"))
    
    # Predict on the test data and add residual noise for stochasticity
    predictions <- predict(model, newdata = data_test_missing, type = "response")
    stochastic_imputations <- rpois(length(predictions), lambda = predictions)
    
    # Ensure imputations are within the original variable range
    min_val <- min(data_train_nonmissing[[target_var]], na.rm = TRUE)
    max_val <- max(data_train_nonmissing[[target_var]], na.rm = TRUE)
    stochastic_imputations <- pmin(pmax(stochastic_imputations, min_val), max_val)
    
  } else if (all(data_train[[target_var]] %in% c(0, 1, NA), na.rm = TRUE)) {
    # For binary data, use logistic regression
    model <- glm(as.formula(paste(target_var, "~", paste(predictors, collapse = " + "))),
                 data = data_train_nonmissing, family = binomial(link = "logit"))
    
    # Predict on the test data and add stochastic noise for stochasticity
    predictions <- predict(model, newdata = data_test_missing, type = "response")
    stochastic_imputations <- rbinom(length(predictions), size = 1, prob = predictions)
    
  } else {
    # For continuous data, use linear regression
    model <- lm(as.formula(paste(target_var, "~", paste(predictors, collapse = " + "))),
                data = data_train_nonmissing)
    
    # Predict on the test data and add residual noise for stochasticity
    residual_sd <- sd(residuals(model))
    predictions <- predict(model, newdata = data_test_missing)
    stochastic_imputations <- predictions + rnorm(length(predictions), mean = 0, sd = residual_sd)
  }
  
  # Replace missing values in the test data with stochastic imputations
  data_test[[target_var]][is.na(data_test[[target_var]])] <- stochastic_imputations
  
  return(data_test)
}


predictors_without_missing <- c("PAT_age", "MED_7d_TPN", "MS_is_ICU_unit", 
                                "MED_L2_7d_L01_ANTINEOPLASTIC_AGENTS", 
                                "MED_L2_7d_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE", 
                                "CAT_catheter_type_binary_all_CVC", 
                                "CAT_catheter_type_binary_all_Tunneled_CVC", 
                                "CAT_catheter_type_binary_all_Dialysis_CVC", 
                                "CAT_catheter_type_binary_all_Port_a_cath", 
                                "CAT_catheter_type_binary_all_PICC")


# mixed-effect model approach
impute_missing_values <- function(data_train, data_test, target_var, predictors, group_var, count_vars, binary_vars, linear_vars) {
  # Load required package
  library(lme4)
  library(glmmTMB)

  # Construct the formula dynamically
  formula <- as.formula(paste(target_var, "~", paste(predictors, collapse = " + "), "+ (1 |", group_var, ")"))
    
  # Initialize the model variable
  model <- NULL
  
  # Set control parameters to suppress convergence warnings
  control_options <- glmerControl(check.conv.grad = "ignore", check.conv.singular = "ignore", check.conv.hess = "ignore")
  
  # Choose the appropriate mixed-effects model based on the target variable type
  if (target_var %in% count_vars) {
    # Poisson mixed-effects model (for count data)
    model <- glmer(formula, data = data_train, na.action = na.exclude, family = poisson(link = "log"), control = control_options)
    # Get the min and max of the target variable in data_train
    min_val <- min(data_train[[target_var]], na.rm = TRUE)
    max_val <- max(data_train[[target_var]], na.rm = TRUE)
  } else if (target_var %in% binary_vars) {
    # Generalized mixed-effects model (for binary data)
    # model <- glmer(formula, data = data_train, na.action = na.exclude, family = binomial(link = "logit"), control = control_options)
    model <- glmmTMB(formula, data = data_train, na.action = na.exclude, family = binomial)
  } else if (target_var %in% linear_vars) {
    # Linear mixed-effects model (for continuous data)
    model <- lmer(formula, data = data_train, na.action = na.exclude)
  } else {
    warning(paste("Skipping", target_var, "- not specified in either generalized_vars or linear_vars."))
    next
  }
    
  # Identify rows with missing values in the target variable
  missing_rows <- is.na(data_test[[target_var]])
    
  # Predict missing values using the fitted model
  predictions <- predict(model, newdata = data_test[missing_rows, ], allow.new.levels = TRUE)
    
  # Round predictions if the target variable is in the generalized list
  if (target_var %in% count_vars) {
    predictions <- exp(predictions)
    predictions <- round(predictions)
    predictions <- pmin(pmax(predictions, min_val), max_val)
  } else if (target_var %in% binary_vars) {
    predictions <- 1 / (1 + exp(-predictions))
    predictions <- round(predictions) 
  }
    
  # Replace missing values with predictions
  data_test[[target_var]][missing_rows] <- predictions
  
  return(data_test)
}


# Define which variables are generalized (binary or count) and which are linear (continuous)
count_vars <- c("CAT_lumens_CVC", "CAT_lumens_Tunneled_CVC", "CAT_lumens_PICC", "CARE_NEU_GCS_score_last")
binary_vars <- c("LAB_is_neutropenia", 
                 # "ADM_admission_source_binary_all_Home", 
                 "MS_medical_specialty_bin_drop_Cardiac", 
                 "MS_medical_specialty_bin_drop_Traumatology",
                 "MS_medical_specialty_bin_drop_Pediatrics")
linear_vars <- c("CARE_VS_temperature_max", 
                 "LAB_urea_last_log", 
                 "LAB_creatinine_last_log", 
                 "LAB_bilirubin_last_log", 
                 "LAB_ciclosporin_last_log", 
                 "LAB_pH_last")


# for mice, create LM^2 and interaction terms
# Define a function to create a new variable -- mice
create_new_mice <- function(data, cat_cols, dropped_levels) {
  # data$LAB_urea_last_log <- log(data$LAB_urea_last)  
  # data$LAB_creatinine_last_log <- log(data$LAB_creatinine_last)
  # data$LAB_bilirubin_last_log <- log(data$LAB_bilirubin_last)
  # data$LAB_ciclosporin_last_log <- log(data$LAB_ciclosporin_last)
  
  # transformations of variables
  data$LM_1 <- data$LM/30
  data$LM_2 <- (data$LM/30)^2
  data$MS_is_ICU_unit_1 <-  data$MS_is_ICU_unit * data$LM_1
  data$MS_is_ICU_unit_2 <-  data$MS_is_ICU_unit * data$LM_2
  
  data <- data %>%
    mutate(CAT_lumens_Total = as.numeric(as.character(CAT_lumens_CVC)) + 
             as.numeric(as.character(CAT_lumens_Tunneled_CVC)) +
             as.numeric(as.character(CAT_lumens_Dialysis_CVC)) +
             as.numeric(as.character(CAT_lumens_Port_a_cath)) +
             as.numeric(as.character(CAT_lumens_PICC))) 
  
  data <- data %>%
    mutate(CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath = ifelse(CAT_catheter_type_binary_all_Tunneled_CVC==1|CAT_catheter_type_binary_all_Port_a_cath==1, 1, 0))
  
  data <- data %>% 
    mutate(LAB_is_neutropenia = as.numeric(LAB_WBC_count_last * 1000 < 500 |
                                             LAB_WBC_Neutrophils_last * 1000 < 500))
  data$ADM_admission_source_binary_all_Home <- as.numeric(as.character(data$ADM_admission_source_binary_all_Home))
  data$CARE_NEU_GCS_score_last <- as.numeric(as.character(data$CARE_NEU_GCS_score_last))
  
  # transform categorical column to binary using `make_col_binary_drop`
  bin_model <- data %>%
    make_col_binary_drop(cat_cols, dropped_levels = dropped_levels)
  data_bin <- bin_model$data
  
  # processing the outcome
  data_bin <- data_bin %>%
    mutate(
      type = if_else(
        (type == "CLABSI" | type == "Death" | type == "Discharge") & eventtime > LM + 7,
        "Censored", type
      ),
      id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
      Tstart = LM,
      eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)
    ) %>%
    filter(eventtime > LM) %>% 
    mutate(
      type = case_when(
        type == "Censored" ~ 0,
        type == "CLABSI" ~ 1,
        type == "Death" ~ 2,
        TRUE ~ 3
      )
    )
  
  return(data_bin)
}


######
create_new_mice2 <- function(data, cat_cols, dropped_levels) {
  # data$LAB_urea_last_log <- log(data$LAB_urea_last)  
  # data$LAB_creatinine_last_log <- log(data$LAB_creatinine_last)
  # data$LAB_bilirubin_last_log <- log(data$LAB_bilirubin_last)
  # data$LAB_ciclosporin_last_log <- log(data$LAB_ciclosporin_last)
  
  # transformations of variables
  data$LM_1 <- data$LM/30
  data$LM_2 <- (data$LM/30)^2
  data$MS_is_ICU_unit_1 <-  data$MS_is_ICU_unit * data$LM_1
  data$MS_is_ICU_unit_2 <-  data$MS_is_ICU_unit * data$LM_2
  
  data <- data %>%
    mutate(CAT_lumens_Total = CAT_lumens_CVC + 
             CAT_lumens_Tunneled_CVC +
             CAT_lumens_Dialysis_CVC +
             CAT_lumens_Port_a_cath +
             CAT_lumens_PICC) 
  
  data <- data %>%
    mutate(CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath = ifelse(CAT_catheter_type_binary_all_Tunneled_CVC==1|CAT_catheter_type_binary_all_Port_a_cath==1, 1, 0))
  
  data <- data %>% 
    mutate(LAB_is_neutropenia = as.numeric(LAB_WBC_count_last * 1000 < 500 |
                                             LAB_WBC_Neutrophils_last * 1000 < 500)) 
  data$ADM_admission_source_binary_all_Home <- as.numeric(as.character(data$ADM_admission_source_binary_all_Home))

  # transform categorical column to binary using `make_col_binary_drop`
  bin_model <- data %>%
    make_col_binary_drop(cat_cols, dropped_levels = dropped_levels)
  data_bin <- bin_model$data
  
  # processing the outcome
  data_bin <- data_bin %>%
    mutate(
      type = if_else(
        (type == "CLABSI" | type == "Death" | type == "Discharge") & eventtime > LM + 7,
        "Censored", type
      ),
      id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
      Tstart = LM,
      eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)
    ) %>%
    filter(eventtime > LM) %>% 
    mutate(
      type = case_when(
        type == "Censored" ~ 0,
        type == "CLABSI" ~ 1,
        type == "Death" ~ 2,
        TRUE ~ 3
      )
    )
  
  return(data_bin)
}


######
create_new_mice3 <- function(data) {
  # data$LAB_urea_last_log <- log(data$LAB_urea_last)  
  # data$LAB_creatinine_last_log <- log(data$LAB_creatinine_last)
  # data$LAB_bilirubin_last_log <- log(data$LAB_bilirubin_last)
  # data$LAB_ciclosporin_last_log <- log(data$LAB_ciclosporin_last)
  
  # transformations of variables
  data$LM_1 <- data$LM/30
  data$LM_2 <- (data$LM/30)^2
  data$MS_is_ICU_unit_1 <-  data$MS_is_ICU_unit * data$LM_1
  data$MS_is_ICU_unit_2 <-  data$MS_is_ICU_unit * data$LM_2
  
  data <- data %>%
    mutate(CAT_lumens_Total = CAT_lumens_CVC + 
             CAT_lumens_Tunneled_CVC +
             CAT_lumens_Dialysis_CVC +
             CAT_lumens_Port_a_cath +
             CAT_lumens_PICC) 
  
  data <- data %>%
    mutate(CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath = ifelse(CAT_catheter_type_binary_all_Tunneled_CVC==1|CAT_catheter_type_binary_all_Port_a_cath==1, 1, 0))
  
  data <- data %>% 
    mutate(LAB_is_neutropenia = as.numeric(LAB_WBC_count_last * 1000 < 500 |
                                             LAB_WBC_Neutrophils_last * 1000 < 500)) 
  data$ADM_admission_source_binary_all_Home <- as.numeric(as.character(data$ADM_admission_source_binary_all_Home))
  

  # processing the outcome
  data_bin <- data %>%
    mutate(
      type = if_else(
        (type == "CLABSI" | type == "Death" | type == "Discharge") & eventtime > LM + 7,
        "Censored", type
      ),
      id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
      Tstart = LM,
      eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)
    ) %>%
    filter(eventtime > LM) %>% 
    mutate(
      type = case_when(
        type == "Censored" ~ 0,
        type == "CLABSI" ~ 1,
        type == "Death" ~ 2,
        TRUE ~ 3
      )
    )
  
  return(data_bin)
}



create_new_mice_yx <- function(data, data_y) {
  # data$LAB_urea_last_log <- log(data$LAB_urea_last)  # Example: Generating random values
  # data$LAB_creatinine_last_log <- log(data$LAB_creatinine_last)
  # data$LAB_bilirubin_last_log <- log(data$LAB_bilirubin_last)
  # data$LAB_ciclosporin_last_log <- log(data$LAB_ciclosporin_last)
  
  data$LM_1 <- data$LM/30
  data$LM_2 <- (data$LM/30)^2
  data$MS_is_ICU_unit_1 <-  data$MS_is_ICU_unit * data$LM_1
  data$MS_is_ICU_unit_2 <-  data$MS_is_ICU_unit * data$LM_2
  
  data <- data %>%
    mutate(CAT_lumens_Total = CAT_lumens_CVC + CAT_lumens_Tunneled_CVC + CAT_lumens_Dialysis_CVC + CAT_lumens_Port_a_cath + CAT_lumens_PICC) 
  
  data <- data %>%
    mutate(CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath = ifelse(CAT_catheter_type_binary_all_Tunneled_CVC==1|CAT_catheter_type_binary_all_Port_a_cath==1, 1, 0))
  
  if (!is.logical(data_y)) {
    data$type <- data_y$type
    data$eventtime <- data_y$eventtime
  }
  
  return(data)
}

######
create_new_mice_yx2 <- function(data, data_y) {
  # data$LAB_urea_last_log <- log(data$LAB_urea_last)  
  # data$LAB_creatinine_last_log <- log(data$LAB_creatinine_last)
  # data$LAB_bilirubin_last_log <- log(data$LAB_bilirubin_last)
  # data$LAB_ciclosporin_last_log <- log(data$LAB_ciclosporin_last)
  
  # transformations of variables
  data$LM_1 <- data$LM/30
  data$LM_2 <- (data$LM/30)^2
  data$MS_is_ICU_unit_1 <-  data$MS_is_ICU_unit * data$LM_1
  data$MS_is_ICU_unit_2 <-  data$MS_is_ICU_unit * data$LM_2
  
  data <- data %>%
    mutate(CAT_lumens_Total = CAT_lumens_CVC + 
             CAT_lumens_Tunneled_CVC +
             CAT_lumens_Dialysis_CVC +
             CAT_lumens_Port_a_cath +
             CAT_lumens_PICC) 
  
  data <- data %>%
    mutate(CAT_catheter_type_binary_all_Tunneled_CVC_Port_a_cath = ifelse(CAT_catheter_type_binary_all_Tunneled_CVC==1|CAT_catheter_type_binary_all_Port_a_cath==1, 1, 0))
  
  data <- data %>% 
    mutate(LAB_is_neutropenia = as.numeric(LAB_WBC_count_last * 1000 < 500 |
                                             LAB_WBC_Neutrophils_last * 1000 < 500)) 
  data$ADM_admission_source_binary_all_Home <- as.numeric(as.character(data$ADM_admission_source_binary_all_Home))
  
  if (!is.logical(data_y)) {
    data <- data %>% select(-c("type", "eventtime")) %>%
      left_join(data_y, by = c("functioneelDossierNr", "CAT_catheter_episode", "LM"))
  }
  
  # processing the outcome
  data <- data %>%
    mutate(
      type = if_else(
        (type == "CLABSI" | type == "Death" | type == "Discharge") & eventtime > LM + 7,
        "Censored", type
      ),
      id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_"),
      Tstart = LM,
      eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7)
    ) %>%
    filter(eventtime > LM) %>% 
    mutate(
      type = case_when(
        type == "Censored" ~ 0,
        type == "CLABSI" ~ 1,
        type == "Death" ~ 2,
        TRUE ~ 3
      )
    )
  
  return(data)
}



# for mice, create LM^2 and interaction terms
# variables used for imputation in mice
vars_mice <- c("PAT_age",
               "ADM_admission_source_binary_all_Home",
               "CAT_catheter_type_binary_all_CVC",
               "CAT_catheter_type_binary_all_Tunneled_CVC",
               "CAT_catheter_type_binary_all_Port_a_cath",
               "CAT_catheter_type_binary_all_PICC",
               "CAT_lumens_CVC",
               "CAT_lumens_Tunneled_CVC",
               "CAT_lumens_Dialysis_CVC",
               "CAT_lumens_Port_a_cath",
               "CAT_lumens_PICC",
               "MED_7d_TPN",
               "MED_L2_7d_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE",
               "MED_L2_7d_L01_ANTINEOPLASTIC_AGENTS",
               "CARE_NEU_GCS_score_last",
               "MS_medical_specialty",
               # "MS_medical_specialty_bin_drop_Cardiac",
               # "MS_medical_specialty_bin_drop_Traumatology",
               # "MS_medical_specialty_bin_drop_Pediatrics",
               "MS_is_ICU_unit",
               "CARE_VS_temperature_max",
               "LAB_is_neutropenia",
               "LAB_urea_last_log",
               "LAB_creatinine_last_log",
               "LAB_bilirubin_last_log",
               "LAB_ciclosporin_last_log",
               "LAB_pH_last")



vars_mice3 <- c("PAT_age",
               "ADM_admission_source_binary_all_Home",
               "CAT_catheter_type_binary_all_CVC",
               "CAT_catheter_type_binary_all_Tunneled_CVC",
               "CAT_catheter_type_binary_all_Port_a_cath",
               "CAT_catheter_type_binary_all_PICC",
               "CAT_lumens_CVC",
               "CAT_lumens_Tunneled_CVC",
               "CAT_lumens_Dialysis_CVC",
               "CAT_lumens_Port_a_cath",
               "CAT_lumens_PICC",
               "MED_7d_TPN",
               "MED_L2_7d_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE",
               "MED_L2_7d_L01_ANTINEOPLASTIC_AGENTS",
               "CARE_NEU_GCS_score_last",
               # "MS_medical_specialty",
               "MS_medical_specialty_bin_drop_Cardiac",
               "MS_medical_specialty_bin_drop_Traumatology",
               "MS_medical_specialty_bin_drop_Pediatrics",
               "MS_is_ICU_unit",
               "CARE_VS_temperature_max",
               "LAB_is_neutropenia",
               "LAB_urea_last_log",
               "LAB_creatinine_last_log",
               "LAB_bilirubin_last_log",
               "LAB_ciclosporin_last_log",
               "LAB_pH_last")

library(brms)

impute_missing_values_multivariate <- function(data_train, data_test, count_vars, binary_vars, linear_vars, predictors, group_var) {
  # Construct formulas for each variable type
  # bf_list <- list()
  
  linear_bfs <- lapply(linear_vars, function(var) {
    bf(as.formula(paste(var, "~", paste(predictors, collapse = " + "), "+ (1 |", group_var, ")")), family = gaussian())
  })
  
  # Binary variables
  binary_bfs <- lapply(binary_vars, function(var) {
    bf(as.formula(paste(var, "~", paste(predictors, collapse = " + "), "+ (1 |", group_var, ")")), family = bernoulli())
  })
  
  # Count variables
  count_bfs <- lapply(count_vars, function(var) {
    bf(as.formula(paste(var, "~", paste(predictors, collapse = " + "), "+ (1 |", group_var, ")")), family = poisson())
  })
  
  # Combine all formulas into a multivariate model
  multivariate_formula <- Reduce(`+`, c(linear_bfs, binary_bfs, count_bfs))
  
  # Fit the multivariate model
  fit <- brm(
    formula = multivariate_formula,
    data = data_train,
    cores = 4,
    chains = 4,
    iter = 2000,
    control = list(adapt_delta = 0.95)
  )
  
  # Predict missing values and impute in the test set
  for (var in c(linear_vars, binary_vars, count_vars)) {
    # Identify rows with missing values in the target variable
    missing_rows <- is.na(data_test[[var]])
    
    if (any(missing_rows)) {
      # Predict values for missing rows
      predictions <- predict(fit, newdata = data_test[missing_rows, ])
      
      # Transform predictions based on variable type
      if (var %in% binary_vars) {
        # Predictions for binary variables
        probs <- predictions[, "Estimate"] # Predicted probabilities
        imputed_values <- as.integer(probs > 0.5) # Threshold at 0.5
      } else if (var %in% count_vars) {
        # Predictions for count variables
        rates <- exp(predictions[, "Estimate"]) # Poisson rates (inverse of log-link)
        imputed_values <- round(rates) # Round to nearest integer
      } else {
        # Predictions for continuous variables
        imputed_values <- predictions[, "Estimate"]
      }
      
      # Replace missing values
      data_test[[var]][missing_rows] <- imputed_values
    }
  }
  
  return(data_test)
}

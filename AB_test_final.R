# Import packages
library(tidyverse)
library(sandwich)
library(lmtest)
library(plm)  # For fixed effects and random effects models

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Preset variable values
n <- 200
case <- "case_3"
required_significance <- 0.05
dependent_var <- "crd_spent"
use_fe <- FALSE  # Set to TRUE to use Fixed Effects, FALSE to use Random Effects

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Create randomly generated data for demonstration
set.seed(42)

# Generate random data for demonstration
data <- data.frame(
  mbr_id = rep(1:n, each = 2),
  crd_spent = round(rnorm(n * 2, mean = 250, sd = 50), 2),
  during_test = rep(0:1, times = n),
  AB = rep(sample(0:1, n, replace = TRUE), each = 2),
  extra_var1 = rnorm(n * 2, mean = 100, sd = 10),
  extra_var2 = rnorm(n * 2, mean = 50, sd = 5)
)

# Include all variables and deselect the unnecessary ones based on the case
if (case == "case_1") {
  selected_data <- data[, setdiff(colnames(data), c("during_test"))]
  independent_vars <- colnames(selected_data)[!colnames(selected_data) %in% c(dependent_var, "mbr_id", "AB")]
  formula <- as.formula(paste(dependent_var, "~ AB +", paste(independent_vars, collapse = " + ")))
} else if (case == "case_2") {
  selected_data <- data[, setdiff(colnames(data), c("AB"))]
  independent_vars <- colnames(selected_data)[!colnames(selected_data) %in% c(dependent_var, "mbr_id", "during_test")]
  formula <- as.formula(paste(dependent_var, "~ during_test +", paste(independent_vars, collapse = " + ")))
} else if (case == "case_3") {
  selected_data <- data[, colnames(data)]
  independent_vars <- colnames(selected_data)[!colnames(selected_data) %in% c(dependent_var, "mbr_id", "during_test", "AB")]
  formula <- as.formula(paste(dependent_var, "~ during_test * AB +", paste(independent_vars, collapse = " + ")))
} else {
  stop("Invalid case. Please choose 'case_1', 'case_2', or 'case_3'.")
}

# Remove unnecessary df to save memory
rm(data)

# Print the dynamically created formula
print(formula)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### case_1 
if (case == "case_1") {
  # Perform the regression for case_1
  model <- lm(formula, data = selected_data)
  
  # Get robust standard errors and insert them into the model
  # HC3 is used for the reason we have a smaller sample size and 
  # we are particularly concerned about leverage points (extreme value predictors)
  robust_se <- coeftest(model, vcov = vcovHC(model, type = "HC3"))
  
  # Print the results
  print(robust_se)
  
  # Graph the results
  # Extract coefficients and standard errors
  coefficients <- coef(robust_se)
  p_values <- robust_se[, 4]
  std_errors <- sqrt(diag(vcovHC(model, type = "HC3")))
  
  # Calculate predicted values and standard errors for control and treatment groups
  predicted_control <- coefficients[1]  # Intercept
  predicted_treatment <- coefficients[1] + coefficients[2]  # Intercept + AB coefficient
  se_control <- std_errors[1]  # Standard error for intercept
  se_treatment <- sqrt(sum(std_errors^2))  # Standard error for the sum of coefficients
  
  # Create a data frame for the dynamite plot
  dynamite_data <- data.frame(
    Group = c("Control", "Treatment"),
    crd_spent = c(predicted_control, predicted_treatment),
    se = c(se_control, se_treatment)
  )
  
  # Print the table
  print(dynamite_data)
  
  # Create the dynamite plot
  graph <- ggplot(dynamite_data, aes(x = Group, y = crd_spent, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.4, alpha = 0.8) +
    geom_errorbar(aes(ymin = crd_spent - se, ymax = crd_spent + se), width = 0.1, position = position_dodge(0.3)) +
    theme_bw() +
    labs(y = "Average Credit Spent") +
    scale_fill_manual(values = c("lightgrey", "lightgrey")) +
    theme(legend.position = "none")
  
  print(graph)
  
  # Summary for management
  if (p_values[2] < required_significance) {
    message <- paste(
      "The analysis indicates a significant effect of the treatment on the dependent variable.",
      "The treatment increases the average credit spent by", round(coefficients[2], 2), "units."
    )
  } else {
    message <- "The analysis indicates no significant effect of the treatment on the dependent variable."
  }
  
  # Print the summary message
  cat(message)
  
  # Remove unnecessary variables
  rm(model, robust_se, coefficients, p_values, std_errors, 
     predicted_control, predicted_treatment, 
     se_control, se_treatment, dynamite_data, graph, message)
}

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### case_2
if (case == "case_2") {
  if (use_fe) {
    # Perform the fixed effects regression for case_2
    model <- plm(formula, data = selected_data, index = "mbr_id", model = "within")
  } else {
    # Perform the random effects regression for case_2
    model <- plm(formula, data = selected_data, index = "mbr_id", model = "random")
  }
  
  # Get robust standard errors and insert them into the model
  # HC3 is used for the reason we have a smaller sample size and 
  # we are particularly concerned about leverage points (extreme value predictors)
  robust_se <- coeftest(model, vcov = vcovHC(model, type = "HC3"))
  
  # Print the results
  print(robust_se)
  
  # Extract p-values from the robust standard errors
  p_values <- robust_se[, "Pr(>|t|)"]
  
  # Graph the results
  # Extract coefficients and standard errors
  coefficients <- robust_se[, "Estimate"]
  std_errors <- robust_se[, "Std. Error"]
  
  # Calculate predicted values and standard errors for before and after groups
  predicted_before <- coefficients["(Intercept)"]  # Intercept
  predicted_after <- coefficients["(Intercept)"] + coefficients["during_test"]  # Intercept + during_test coefficient
  se_before <- std_errors["(Intercept)"]  # Standard error for intercept
  
  # Extract the variance-covariance matrix
  vcov_matrix <- vcovHC(model, type = "HC3")
  
  # Calculate the standard error for the predicted value after the treatment
  se_after <- sqrt(
    vcov_matrix["(Intercept)", "(Intercept)"] + 
      vcov_matrix["during_test", "during_test"] + 
      2 * vcov_matrix["(Intercept)", "during_test"]
  )
  
  # Create a data frame for the dynamite plot
  dynamite_data <- data.frame(
    Group = c("Before", "After"),
    crd_spent = c(predicted_before, predicted_after),
    se = c(se_before, se_after)
  )
  
  # Print the table
  print(dynamite_data)
  
  # Create the dynamite plot
  graph <- ggplot(dynamite_data, aes(x = Group, y = crd_spent, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.4, alpha = 0.8) +
    geom_errorbar(aes(ymin = crd_spent - se, ymax = crd_spent + se), width = 0.1, position = position_dodge(0.3)) +
    theme_bw() +
    labs(y = "Average Credit Spent") +
    scale_fill_manual(values = c("lightgrey", "lightgrey")) +
    theme(legend.position = "none")
  
  print(graph)
  
  # Summary for management
  if (p_values["during_test"] < required_significance) {
    message <- paste(
      "The analysis indicates a significant effect of the treatment on the dependent variable.",
      "The treatment increases the average credit spent by", round(coefficients["during_test"], 2), "units."
    )
  } else {
    message <- "The analysis indicates no significant effect of the treatment on the dependent variable."
  }
  
  # Print the summary message
  cat(message)
  
  # Remove unnecessary variables
  rm(coefficients, std_errors, predicted_before, predicted_after, 
     se_before, vcov_matrix, se_after, dynamite_data, graph, message)
}

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### case_3
if (case == "case_3") {
  if (use_fe) {
    # Perform the fixed effects regression for case_3
    model <- plm(formula, data = selected_data, index = "mbr_id", model = "within")
  } else {
    # Perform the random effects regression for case_3
    model <- plm(formula, data = selected_data, index = "mbr_id", model = "random")
  }
  
  # Get robust standard errors and insert them into the model
  # HC3 is used for the reason we have a smaller sample size and 
  # we are particularly concerned about leverage points (extreme value predictors)
  robust_se <- coeftest(model, vcov = vcovHC(model, type = "HC3"))
  
  # Print the results
  print(robust_se)
  
  # Extract p-values from the robust standard errors
  p_values <- robust_se[, "Pr(>|t|)"]
  
  # Extract coefficients
  coefficients <- robust_se[, "Estimate"]
  
  # Parameter-based calculations
  beta_0 <- coefficients["(Intercept)"]
  beta_1 <- coefficients["AB"]
  beta_2 <- coefficients["during_test"]
  beta_3 <- coefficients["during_test:AB"]
  
  control_before <- beta_0
  control_after <- beta_0 + beta_2
  treated_before <- beta_0 + beta_1
  treated_after <- beta_0 + beta_1 + beta_2 + beta_3
  
  # Create a data frame for the parameter-based calculations table
  param_table <- data.frame(
    Condition = c("Before", "After"),
    Control = c(control_before, control_after),
    Treated = c(treated_before, treated_after)
  )
  
  # Print the parameter-based calculations table
  print(param_table)
  
  # Summary for management
  if (p_values["during_test:AB"] < required_significance) {
    message <- paste(
      "The analysis indicates a significant interaction effect of the treatment and time on the dependent variable.",
      "The interaction effect increases the average credit spent by", round(coefficients["during_test:AB"], 2), "units."
    )
  } else {
    message <- "The analysis indicates no significant interaction effect of the treatment and time on the dependent variable."
  }
  
  # Print the summary message
  cat(message)
  
  # Remove created variables to clean up environment
  rm(coefficients, beta_0, beta_1, beta_2, beta_3, 
     control_before, control_after, treated_before, treated_after, 
     param_table, message)
}

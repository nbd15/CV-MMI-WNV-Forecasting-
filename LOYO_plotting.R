# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(data.table)
  library(rio)
  library(caret)
})

source("glmm_functions.R")

n_preds <- 4
figure_folder <- "figures/loyo_plots/"
nldas_mle_annual <- read.csv("data/cv_n500_mle_2006_2022_over_11_obs.csv")

# get_ens_weights
#
# Inputs:
#   1. start_year, end_year -- years used to get LOYO model results
#   2. year -- any year between start_year and end_year
#   3. month_range -- which months were used to define the LOYO models
#
# Outputs:
#   1. model_wts -- adjusted weights of models in ensemble
get_ens_weights <- function(start_year, end_year, year, month_range) {
  
  fname = paste0("model_outputs/2006-2021/cv_", start_year, "_", end_year, "_ensemble_results_", year, ".xlsx")
  
  data_lst <- import_list(fname)
  model_ests <- setDT(data_lst$`Model Parameters`)
  
  model_wts <- get_model_weights(model_ests, month_range = month_range)
  
  if (is.character(model_wts)) {
    
    filler_output <- data.frame(matrix(c(0, 0, 0, 0, 0, 0, "tmp_1", 0, year), nrow = 1))
    
    colnames(filler_output) <- c("model_num", "adj.wts", "Estimate", "Std. Error", 
                                 "z value", "pval", "varnames", "aicc", "year")
    
    filler_output$adj.wts <- as.numeric(filler_output$adj.wts)
    filler_output$year <- as.numeric(filler_output$year)
    
    return(filler_output)
    
  }
  
  model_wts <- model_wts %>%
    mutate(year = year)
  
  return(model_wts)
  
}

start_year <- 2006
end_year <- 2021

# get ensemble weights for all models between 2006 and 2021
ens_weights <- rbindlist(lapply(start_year:end_year, get_ens_weights, 
                                start_year = start_year, end_year = end_year, month_range = c(11, 7)))

# calculate the weighting of each environmental variable within an ensemble
env_var_weights <- ens_weights %>%
  filter(varnames != "(Intercept)") %>%
  mutate(temp_wts = adj.wts / 4) %>%
  group_by(year, varnames) %>%
  summarize(plot_wts = sum(temp_wts)) %>%
  separate(col = varnames, into = c("variable", "month"), sep = "_") %>%
  mutate(variable = ifelse(variable == "tmp", "Temp", "ET"),
         month = month.abb[as.numeric(month)],
         plot_varname = paste(month, variable),
         plot_varname = factor(plot_varname, levels = c(paste(month.abb[c(11:12, 1:7)], "ET"),
                                                        paste(month.abb[c(11:12, 1:7)], "Temp"))))

# plot variable weighting by year
var_weights_plot <- env_var_weights %>%
  ggplot(aes(x = year, y = plot_wts, fill = plot_varname)) +
  geom_col() + 
  labs(x = "Year",
       y = "Variable Importance",
       fill = "Environmental\nVariables",
       title = paste0("LOYO Variable Weights: ", start_year, " to ", end_year, "\n")) + 
  scale_fill_manual(limits = levels(env_var_weights$plot_varname),
                    values = c(colorRampPalette(c("#03018C", "#9EC2FF"))(9),
                               colorRampPalette(c("#DC1C13", "#F6BDC0"))(9))) + 
  guides(fill = guide_legend(ncol = 2, byrow = F)) + 
  scale_x_continuous(name = "Year", breaks = seq(start_year, end_year, 2))

var_weights_plot

ggsave(paste0(figure_folder, "var_weights_", start_year, "_", end_year, ".png"),
       var_weights_plot, width = 12, height = 6)

################################################################################

# Sensitivity and Specificity

# get_ens_ests
#
# Inputs:
#   1. start_year, end_year -- years used to get LOYO model results
#   2. year -- any year between start_year and end_year
#   3. threshold -- threshold for determining agreement between observed and predicted
#
# Outputs:
#   1. a dataframe containing the observed infection rate, ensemble estimate,
#      best model estimate, and booleans showing if they are greater than the
#      threshold
get_ens_ests <- function(start_year, end_year, year, threshold) {
  
  # import LOYO model results
  fname <- paste0("model_outputs/2006-2021/cv_", start_year, "_", end_year, "_ensemble_results_", year, ".xlsx")
  
  data_lst <- import_list(fname)
  
  # extract ensemble and best model estimates, select columns
  ens_ests <- data_lst$`Ensemble Estimates`
  best_ests <- data_lst$`Best Model Estimates`
  
  ens_ests <- ens_ests %>%
    dplyr::select(nldasID, ensemble_est)
  
  best_ests <- best_ests %>%
    dplyr::select(nldasID, est)
  
  # combine ensemble and best model estimates
  output_ests <- left_join(ens_ests, best_ests, by = "nldasID") %>%
    mutate(year = year) %>%
    rename("best_est" = est)
  
  # combine estimates with observed
  output_ests <- left_join(output_ests, nldas_mle_annual, by = c("nldasID", "year")) %>%
    dplyr::select(nldasID:year, MLE) %>%
    rename("year_observed" = year,
           "observed" = MLE)
  
  # create boolean variables for whether value is greater than or equal to threshold
  output_ests <- output_ests %>%
    mutate(obs_gt_thres = factor(observed >= threshold, levels = c(TRUE, FALSE)),
           best_est_gt_thres = factor(best_est >= threshold, levels = c(TRUE, FALSE)),
           ens_est_gt_thres = factor(ensemble_est >= threshold, levels = c(TRUE, FALSE)))
  
}

# get_sens_spec
#
# Inputs:
#   1. start_year, end_year -- years used to get LOYO model results
#   2. threshold -- threshold for determining agreement between observed and predicted
#
# Outputs:
#   1. a plot of observed v. ensemble and best model estimates
#   2. a dataframe containing the sensitivity, specificity, and RMSE of the 
#      ensemble and best models
get_sens_spec <- function(start_year, end_year, threshold) {
  
  # get combined ensemble estimates for each year between start_year and end_year
  predictions_observed <- rbindlist(lapply(start_year:end_year, get_ens_ests, start_year = start_year, end_year = end_year, threshold = threshold))
  write.csv(predictions_observed, paste0(figure_folder, "obs_pred_loyo_", start_year, "_", end_year, "_", threshold, "_thres.csv"), row.names = FALSE)
  
  # set colors for plotting
  colors <- c("Observed" = "#06ABEB", "Ensemble Estimate" = "orange", "Best Model Estimate" = "#DC298D")
  
  # sort by observed infection rate, plot alongside ensemble and best estimates
  loyo_plot <- predictions_observed %>%
    filter(!is.na(observed)) %>%
    arrange(observed) %>%
    mutate(rnum = 1:nrow(predictions_observed %>% filter(!is.na(observed)))) %>%
    ggplot() + 
    geom_line(aes(x = rnum, y = observed, color = "Observed")) + 
    geom_point(aes(x = rnum, y = ensemble_est, color = "Ensemble Estimate"), shape = 2) + 
    geom_point(aes(x = rnum, y = best_est, color = "Best Model Estimate"), shape = 4) +
    geom_hline(yintercept = threshold, color = "#1f7d1f") + 
    labs(x = expression("Observations Sorted By" ~ I[M]), 
         y = expression(I[M])) +
    theme(legend.title = element_blank()) +
    scale_color_manual(values = colors)
  
  ggsave(paste0(figure_folder, "sens_spec_", start_year, "_", end_year, "_", threshold, "_thres.png"), loyo_plot,
         width = 6, height = 6)
  
  # calculate sensitivity and specificity of observed v. best model estimates
  obs_v_best <- confusionMatrix(predictions_observed$obs_gt_thres, predictions_observed$best_est_gt_thres)
  obs_v_best_sens <- round(obs_v_best$byClass['Sensitivity'], 4)
  obs_v_best_spec <- round(obs_v_best$byClass['Specificity'], 4)
  
  # calculate sensitivity and specificity of observed v. ensemble estimates
  obs_v_ens <- confusionMatrix(predictions_observed$obs_gt_thres, predictions_observed$ens_est_gt_thres)
  obs_v_ens_sens <- round(obs_v_ens$byClass['Sensitivity'], 4)
  obs_v_ens_spec <- round(obs_v_ens$byClass['Specificity'], 4)
  
  # calculate RMSE of observed v. estimates
  obs_v_best_rmse <- sqrt(mean((predictions_observed$observed - predictions_observed$best_est)^2, na.rm = TRUE))
  obs_v_ens_rmse <- sqrt(mean((predictions_observed$observed - predictions_observed$ensemble_est)^2, na.rm = TRUE))
  
  df = data.frame(model_type = "LOYO", start_year, end_year, obs_v_best_sens, 
                  obs_v_best_spec, obs_v_best_rmse, obs_v_ens_sens, obs_v_ens_spec, obs_v_ens_rmse)
  
  return(df)
  
}

all_models_sens_spec <- get_sens_spec(start_year = 2006, end_year = 2021, threshold = 1)

write.csv(all_models_sens_spec, paste0(figure_folder, "sensitivity_specificity.csv"), row.names = FALSE)

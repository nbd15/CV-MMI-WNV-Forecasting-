# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(data.table)
  library(rio)
  library(caret)
})

# Used: get_model_weights, prediction_from_model_wilcox
source("glmm_functions.R")

n_preds <- 4
figure_folder <- "figures/loyo_plots/"
nldas_mle_annual <- read.csv("data/cv_n500_mle_2006_2022_over_11_obs.csv")

start_year <- 2006
end_year <- 2021

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
  
  fname = paste0("model_outputs/", start_year, "-", end_year, "/cv_", 
                 start_year, "_", end_year, "_ensemble_results_", year, ".xlsx")
  
  data_lst <- import_list(fname)
  model_ests <- setDT(data_lst$`Model Parameters`)
  
  model_wts <- get_model_weights(model_ests, month_range = month_range)
  
  model_wts <- model_wts %>%
    mutate(year = year)
  
  return(model_wts)
  
}

# get ensemble weights for all models between 2006 and 2021
ens_weights <- rbindlist(lapply(start_year:end_year, get_ens_weights, 
                                start_year = start_year, end_year = end_year, month_range = c(10, 7)))

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
         plot_varname = factor(plot_varname, levels = c(paste(month.abb[c(10:12, 1:7)], "ET"),
                                                        paste(month.abb[c(10:12, 1:7)], "Temp")))) %>%
  ungroup() %>%
  group_by(year) %>%
  arrange(year, plot_varname) %>%
  mutate(cum_wt = cumsum(plot_wts),
         text_y = 1 - (cum_wt - plot_wts / 2)) %>%
  ungroup() %>%
  group_by(plot_varname) %>%
  mutate(first_obs = row_number(plot_varname) == 1) %>%
  ungroup() %>%
  mutate(bar_text = ifelse(first_obs, as.character(plot_varname), ""))

# plot variable weighting by year
var_weights_plot <- env_var_weights %>%
  ggplot(aes(x = year, y = plot_wts, fill = plot_varname)) +
  geom_col() + 
  labs(x = "Year",
       y = "Variable Importance",
       #title = paste0("LOYO Variable Weights: ", start_year, " to ", end_year, "\n"),
       fill = "Environmental\nVariables") + 
  scale_fill_manual(limits = levels(env_var_weights$plot_varname),
                    values = c(colorRampPalette(c("#03018C", "#9EC2FF"))(10),
                               colorRampPalette(c("#DC1C13", "#F6BDC0"))(10))) + 
  guides(fill = guide_legend(ncol = 2, byrow = F)) + 
  scale_x_continuous(name = "Year", breaks = seq(start_year, end_year, 2)) +
  theme_minimal() + 
  geom_text(aes(x = year, y = text_y, label = bar_text), size = 4, color = "white") + 
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

var_weights_plot

ggsave(paste0(figure_folder, "new_var_weights_", start_year, "_", end_year, ".jpg"),
       var_weights_plot, width = 16, height = 11)

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
  
  fname = paste0("model_outputs/", start_year, "-", end_year, "/cv_", 
                 start_year, "_", end_year, "_ensemble_results_", year, ".xlsx")
  
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


get_null_mle <- function(i_year) {
  
  null_mle <- nldas_mle_annual %>%
    filter(year %in% start_year:end_year & year != i_year) %>%
    group_by(nldasID) %>%
    summarize(null_mle_val = mean(MLE)) %>%
    mutate(year_observed = i_year)
  
  return(null_mle)
  
}

null_mle_all_years <- do.call("rbind", lapply(start_year:end_year, get_null_mle))

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
  predictions_observed <- left_join(predictions_observed, null_mle_all_years)
  
  predictions_observed <- predictions_observed %>%
    mutate(null_gt_thres = factor(null_mle_val >= threshold, levels = c(TRUE, FALSE))) %>%
    filter(!is.na(observed))
  
  write.csv(predictions_observed, paste0(figure_folder, "obs_pred_loyo_", start_year, "_", end_year, "_", threshold, "_thres.csv"), row.names = FALSE)
  
  # set colors for plotting
  colors <- c("Observed" = "#06ABEB", "Ensemble Estimate" = "orange", "Best Model Estimate" = "#DC298D", "Null Model Estimate" = "black")
  
  # sort by observed infection rate, plot alongside ensemble and best estimates
  loyo_plot <- predictions_observed %>%
    filter(!is.na(observed)) %>%
    arrange(observed, ensemble_est) %>%
    mutate(rnum = 1:nrow(predictions_observed %>% filter(!is.na(observed)))) %>%
    ggplot() + 
    geom_line(aes(x = rnum, y = observed, color = "Observed")) + 
    geom_point(aes(x = rnum, y = ensemble_est, color = "Ensemble Estimate"), shape = 2) + 
    geom_point(aes(x = rnum, y = best_est, color = "Best Model Estimate"), shape = 4) +
    geom_point(aes(x = rnum, y = null_mle_val, color = "Null Model Estimate"), shape = 6) + 
    geom_hline(yintercept = threshold, color = "#1f7d1f") + 
    labs(x = expression("Observations Sorted By" ~ I[M]), 
         y = expression(I[M])) +
    scale_color_manual(values = colors) + 
    scale_y_continuous(breaks = c(0, threshold, 5, 10, 15, 20)) + 
    theme_minimal() + 
    theme(legend.title = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(paste0(figure_folder, "sens_spec_", start_year, "_", end_year, "_", threshold, "_thres.jpg"), loyo_plot,
         width = 10, height = 6)
  
  # calculate sensitivity and specificity of observed v. best model estimates
  obs_v_best <- confusionMatrix(reference = predictions_observed$obs_gt_thres, data = predictions_observed$best_est_gt_thres, positive = "TRUE")
  obs_v_best_sens <- round(obs_v_best$byClass['Sensitivity'], 4)
  obs_v_best_spec <- round(obs_v_best$byClass['Specificity'], 4)
  
  # calculate sensitivity and specificity of observed v. ensemble estimates
  obs_v_ens <- confusionMatrix(reference = predictions_observed$obs_gt_thres, data = predictions_observed$ens_est_gt_thres, positive = "TRUE")
  obs_v_ens_sens <- round(obs_v_ens$byClass['Sensitivity'], 4)
  obs_v_ens_spec <- round(obs_v_ens$byClass['Specificity'], 4)
  
  # calculate sensitivity and specificity of observed v. ensemble estimates
  obs_v_null <- confusionMatrix(reference = predictions_observed$obs_gt_thres, data = predictions_observed$null_gt_thres, positive = "TRUE")
  obs_v_null_sens <- round(obs_v_null$byClass['Sensitivity'], 4)
  obs_v_null_spec <- round(obs_v_null$byClass['Specificity'], 4)
  
  
  # calculate RMSE of observed v. estimates
  obs_v_best_rmse <- sqrt(mean((predictions_observed$observed - predictions_observed$best_est)^2, na.rm = TRUE))
  obs_v_ens_rmse <- sqrt(mean((predictions_observed$observed - predictions_observed$ensemble_est)^2, na.rm = TRUE))
  obs_v_null_rmse <- sqrt(mean((predictions_observed$observed - predictions_observed$null_mle_val)^2, na.rm = TRUE))
  
  df = data.frame(model_type = "LOYO", start_year, end_year, 
                  obs_v_best_sens, obs_v_best_spec, obs_v_best_rmse, 
                  obs_v_ens_sens, obs_v_ens_spec, obs_v_ens_rmse,
                  obs_v_null_sens, obs_v_null_spec, obs_v_null_rmse)
  
  return(df)
  
}

all_models_sens_spec <- get_sens_spec(start_year, end_year, threshold = 1)

write.csv(all_models_sens_spec, paste0(figure_folder, start_year, "_", end_year, "_sensitivity_specificity.csv"), row.names = FALSE)

################################################################################

# Wilcox Signed Rank Test

cv_nldas_dt <- fread("data/cv_monthly_nldas_2005_2022.csv")

nldas_cells <- unique(nldas_mle_annual$nldasID)

nldas_dt <- cv_nldas_dt %>%
  filter(nldasID %in% nldas_cells)

get_predictions_for_year_wilcox <- function(i_year) {
  
  fname = paste0("model_outputs/", start_year, "-", end_year, "/cv_", 
                 start_year, "_", end_year, "_ensemble_results_", i_year, ".xlsx")
  
  model_ests <- setDT(read_excel(fname, sheet = "Model Parameters"))
  
  model_wts <- get_model_weights(model_ests, month_range = c(10, 7))
  
  nldas_dt[,nldasID := factor(nldasID)]
  nldas_ids <- unique(nldas_dt$nldasID)
  nldas_dt <- nldas_dt[order(nldasID, year, month)]
  
  # original and new column names for lag
  og_cols <- c("temp", "evpsfc")
  new_cols <- c("tmp_lag", "evp_lag")
  
  # lag environmental variables by 4 months
  # each year will now have data from September to August instead of January to December
  env_dt <- nldas_dt[, (new_cols) := data.table::shift(.SD, 4), by = nldasID, .SDcols = og_cols]
  
  # calculate standardized environmental values
  env_dt_sub <- env_dt[year %in% start_year:end_year & year != i_year]
  env_std <- env_dt_sub[,.(year,tmp_std = scale(tmp_lag, scale = T, center = T)[,1],
                           evp_std = scale(evp_lag, scale = T, center = T)[,1]),.(nldasID,month)]
  
  # calculate mean and sd for environmental variables
  # used in scaling test data
  env_means <- env_dt[year %in% start_year:end_year & year != i_year,.(mean_tmp = mean(tmp_lag,na.rm = T), sd_tmp = sd(tmp_lag,na.rm = T), 
                                                                       mean_evp = mean(evp_lag,na.rm = T), sd_evp = sd(evp_lag,na.rm = T)),.(nldasID, month)]
  
  # pivot standardized environmental vars from long to wide
  env_wide <- dcast(env_std, nldasID + year ~ month, value.var = c("tmp_std","evp_std"))
  setnames(env_wide, c("nldasID","year", paste0("tmp_",c(9:12,1:8)), paste0("evp_",c(9:12,1:8))))
  
  # join wide standardized environmental data to MLE
  train_dt <- env_wide[nldas_mle_annual, on = c("nldasID","year"),nomatch = 0]
  
  train_dt$nldasID <- as.factor(train_dt$nldasID)
  
  # generate all combos of env variables for modeling
  var.dt <- as.data.table(t(combn(names(train_dt)[3:26], n_preds)))
  
  # find models in ensemble
  ens_models_list <- unique(model_wts$model_num)
  
  # create test dataset standardized by the training means and standard deviations
  test_dt <- env_means[env_dt[year == i_year], on = c("nldasID","month"), nomatch = 0]
  test.std <- test_dt[,.(tmp_std = (tmp_lag - mean_tmp)/sd_tmp, evp_std = (evp_lag - mean_evp)/sd_evp),.(nldasID,month,year)]
  test.std[is.na(test.std)] <- 0
  test_dt_w_env_data <- dcast(test.std, nldasID + year ~ month, value.var = c("tmp_std","evp_std"))
  setnames(test_dt_w_env_data, c("nldasID","year", paste0("tmp_",c(9:12,1:8)), paste0("evp_",c(9:12,1:8))))
  
  env_means_names <- c("nldasID", "month", "mean_tmp", "mean_evp")
  
  # create test dataset with no added environmental data
  raw_test_dt <- env_means[, ..env_means_names]
  colnames(raw_test_dt) <- sub("mean", "lag", colnames(raw_test_dt))
  test_dt <- env_means[raw_test_dt, on = c("nldasID", "month")]
  test.std <- test_dt[,.(tmp_std = (lag_tmp - mean_tmp)/sd_tmp, evp_std = (lag_evp - mean_evp)/sd_evp),.(nldasID,month)]
  test_dt_no_env_data <- dcast(test.std, nldasID ~ month, value.var = c("tmp_std","evp_std"))
  setnames(test_dt_no_env_data, c("nldasID", paste0("tmp_",c(9:12,1:8)), paste0("evp_",c(9:12,1:8))))
  test_dt_no_env_data <- test_dt_no_env_data %>%
    mutate(year = i_year)
  
  # get predicted values for models in ensemble
  nldas_4p_all_model_fcst <- rbindlist(lapply(ens_models_list, 
                                              prediction_from_model_wilcox, 
                                              train_dt = train_dt,
                                              test_dt_env = test_dt_w_env_data,
                                              test_dt_no_env = test_dt_no_env_data,
                                              var.dt = var.dt,
                                              m.wts = model_wts), idcol = "model_num")
  
  nldas_4p_all_model_fcst[,model_num := ens_models_list[model_num]]
  
  # combine model estimates using weighted sum, calculate estimate variance
  nldas_4p_ens_model_fcst <- unique(nldas_4p_all_model_fcst[,.(est_env, est_no_env,
                                                               var_env, var_no_env, 
                                                               weight,
                                                               ensemble_est_env = sum(weight*est_env),
                                                               ensemble_est_no_env = sum(weight * est_no_env)),.(nldasID)][
                                                                 ,.(ensemble_est_env,
                                                                    ensemble_est_no_env,
                                                                    ensemble_var_env = (sum(weight*sqrt(var_env + (ensemble_est_env - est_env)^2)))^2,
                                                                    ensemble_var_no_env = (sum(weight*sqrt(var_no_env + (ensemble_est_no_env - est_no_env)^2)))^2),nldasID])
  
  
  # get null model prediction for given year
  null_mle <- nldas_mle_annual %>%
    filter(year != i_year & year %in% start_year:end_year) %>%
    group_by(nldasID) %>%
    summarize(null_mle_pred = mean(MLE))
  
  # get observed infection rate
  obs_year_i <- nldas_mle_annual %>%
    filter(year == i_year) %>%
    dplyr::select(year, nldasID, MLE)
  
  # combine all predictions and calculate absolute error between prediction and observed infection rates
  year_i_preds <- left_join(left_join(obs_year_i, nldas_4p_ens_model_fcst %>% dplyr::select(nldasID, ensemble_est_env, ensemble_est_no_env)), null_mle)
  
  year_i_preds <- year_i_preds %>%
    mutate(ensemble_env_err = abs(MLE - ensemble_est_env),
           ensemble_no_env_err = abs(MLE - ensemble_est_no_env),
           null_mle_err = abs(MLE - null_mle_pred))
  
  return(year_i_preds)
  
}

all_years_pred_errs <- do.call("rbind", lapply(start_year:end_year,
                                               get_predictions_for_year_wilcox))

write.csv(all_years_pred_errs, paste0(figure_folder, "wilcox_errors_", start_year, "_", end_year, ".csv"), row.names = FALSE)

wilcox.test(all_years_pred_errs$ensemble_env_err, all_years_pred_errs$ensemble_no_env_err, paired = TRUE, alternative = "less")
wilcox.test(all_years_pred_errs$ensemble_env_err, all_years_pred_errs$null_mle_err, paired = TRUE, alternative = "less")


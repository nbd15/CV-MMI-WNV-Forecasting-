# glmm_functions.R contains nldas_prediction_function, which uses glmm_wts and admb_fun
source("glmm_functions.R")

# library for Excel output
library(writexl)

#######################################################################################
set.seed(153)

# load data, convert nldasID to factor
nldas_mle_annual <- fread("data/cv_n500_mle_2006_2022_over_11_obs.csv")
nldas_mle_annual[,nldasID := factor(nldasID)]
nldas_dt <- fread("data/cv_monthly_nldas_2005_2022.csv")

# include cells with substantial trapping effort
nldas_dt <- nldas_dt[nldasID %in% unique(nldas_mle_annual$nldasID)]
nldas_dt[,nldasID := factor(nldasID)]
nldas_ids <- unique(nldas_dt$nldasID)

nldas_dt <- nldas_dt[order(nldasID, year, month)]

# range of years and months to include in ensemble
start_year <- 2006
end_year <- 2021
start_month <- 11
end_month <- 7

years <- start_year:end_year

# original and new column names for lag
og_cols <- c("temp", "evpsfc")
new_cols <- c("tmp_lag", "evp_lag")

# shift environmental data by 4 months
# each year now represented as September to August
env_dt <- nldas_dt[, (new_cols) := shift(.SD, 4), by = nldasID, .SDcols = og_cols]
env_dt <- env_dt[year %in% years]

for (i.year in start_year:end_year) {
  
  # standardize training data
  env_std <- env_dt[year != i.year,.(year,
                                     tmp_std = scale(tmp_lag, scale = T, center = T)[,1],
                                     evp_std = scale(evp_lag, scale = T, center = T)[,1]),.(nldasID,month)]
  
  # save means and standard deviations for standardizing test data
  env_means <- env_dt[year != i.year,.(mean_tmp = mean(tmp_lag),
                                       sd_tmp = sd(tmp_lag),
                                       mean_evp = mean(evp_lag),
                                       sd_evp = sd(evp_lag)),.(nldasID, month)]
  
  # pivot from long to wide
  env_wide <- dcast(env_std, nldasID + year ~ month, value.var = c("tmp_std","evp_std"))
  
  setnames(env_wide, c("nldasID","year",
                       paste0("tmp_",c(9:12,1:8)),
                       paste0("evp_",c(9:12,1:8))))
  
  # combine environmental and mosquito infection data
  train_dt <- env_wide[nldas_mle_annual, on = c("nldasID","year"),nomatch = 0]
  
  # Get test data
  test.dt <- env_dt[year == i.year]
  test.dt <- env_means[test.dt, on = c("nldasID","month")]
  
  # standardize test data using training data means
  test.std <- test.dt[,.(tmp_std = (tmp_lag - mean_tmp)/sd_tmp,
                         evp_std = (evp_lag - mean_evp)/sd_evp),.(nldasID,month,year)]
  
  # pivot from long to wide
  test_dt <- dcast(test.std, nldasID + year ~ month, 
                   value.var = c("tmp_std","evp_std"))
  setnames(test_dt, c("nldasID","year",
                      paste0("tmp_",c(9:12,1:8)),
                      paste0("evp_",c(9:12,1:8))))
  n_preds <- 4
  
  # generate ensemble predictions and model parameters
  out_4p <- nldas_prediction_function(train.dt = train_dt,
                                      test.dt = test_dt, n_predictors = n_preds,
                                      varfrom = 3, varto = 26, ncores = 32, 
                                      month_range = c(start_month, end_month))
  
  output <- list("Ensemble Estimates" = out_4p$final.est,
                 "Model Parameters" = out_4p$models,
                 "Model Weights" = out_4p$model_weights,
                 "Best Model Estimates" = out_4p$best.est)
  
  write_xlsx(output, paste0("model_outputs/", start_year, "_", end_year, "/cv_", start_year, "_", end_year,"_ensemble_results_", i.year, ".xlsx", sep = ""))
  
}

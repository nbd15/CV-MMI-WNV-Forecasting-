#######################         Introduction         ###########################

# NLDAS Forecasts and Plotting
#
# Inputs:
#   1. Google Maps API key
#   2. Shapefile for region
#   3. Shapefile for NLDAS grids of region
#   4. Bounding box coordinates for region of interest -- used to generate terrain object
#   5. All model parameters from model output -- for predictions and plotting
#   6. Model input files (MLE and NLDAS)
#   7. Number of predictors in each model
#
# Functions:
#   1. get_model_weights -- generates ensemble model weights for the specified subset of months
#   2. get_forecast_for_year -- generates an ensemble prediction for a given year
#      by compiling the results of individual models in the ensemble
#   3. plot_observed_predicted_agreement -- plots observed and predicted infection rates for 2022
#      and the agreement between observed and predicted
#   4. get_monthly_predictions -- generates a real-time ensemble forecast for months in 2022
#   5. plot_observed_predicted_agreement_monthly -- plots 2022 observed infection rate
#      and real-time monthly predictions with agreement between observed and predicted
#   6. get_rel_imp -- calculates relative importance of each environmental variable in ensemble
#   7. plot_scatterpie -- generates scatterpie plot for Nov-Jul ensemble model subset

# load packages
suppressPackageStartupMessages({
  library(lubridate)
  library(sf)
  library(raster)
  library(fst)
  library(stringr)
  library(parallel)
  library(ggplot2)
  library(ggmap)
  library(gridExtra)
  library(ggspatial)
  library(mapview)
  library(dplyr)
  library(scatterpie)
  library(leaflet)
  library(MMWRweek)
  library(readxl)
  library(tidyverse)
  library(grid)
  library(ggpubr)
  library(data.table)
})

# change global options
theme_set(theme_classic())
sf_use_s2(FALSE)
options(scipen = 999)

# get functions from glmm_functions
# Used: glmm_wts, prediction_from_model, get_model_weights
source("glmm_functions.R")

#######################            Inputs            ###########################

# Input Values

# ggmaps Terrain
# insert Google Maps API key here
# register_google(key = ...)

# number of predictors in each model
# 4 predictors used in provided model output
n_preds <- 4

# set figure output folder name and MLE definition
figure_folder <- "figures/"
MLE_def_str <- "MLE: 500+ mosquitoes tested in at least 12 years"

# import model output parameters
model_ests <- setDT(read_excel("model_outputs/2006-2021/cv_2006_2021_all_years_over_11_obs.xlsx", sheet = "Model Parameters"))
head(model_ests)

# import shapefiles -- all of Coachella Valley and intersecting NLDAS grid cells
region_shp <- st_read("shapefiles/coachella_valley_latlon.shp", quiet = TRUE)
nldas_polygons <- st_read("shapefiles/cv_nldas.shp", quiet = TRUE) %>% rename("nldasID" = NLDAS_ID)

# load NLDAS and MLE data
nldas_dt <- fread("data/cv_monthly_nldas_2005_2022.csv")
nldas_mle_annual <- fread("data/cv_n500_mle_2006_2022_over_11_obs.csv")

# Extract components of shapefile, set terrain
nldas_ids <- nldas_polygons %>% st_drop_geometry() %>% setDT()
terrain <- get_map(location = c(-116.625, 33.375, -115.875, 34.000), zoom = 10, maptype = "terrain")

#######################         Preprocessing        ##########################

# sort by month within each NLDAS cell
nldas_dt <- nldas_dt[order(year,nldasID,month)]
nldas_dt[,nldasID := factor(nldasID)]

# create binned MLE variable
threshold <- round(quantile(nldas_mle_annual$MLE, 0.75)[[1]], 1)
midpoint <- round(threshold / 2, 1)

nldas_mle_annual[,nldasID := factor(nldasID)]
nldas_mle_annual[, mle_cat := cut(MLE, 
                                  breaks = c(0, midpoint, threshold, 1E30),
                                  labels = c(paste0("<", midpoint),
                                             paste0(midpoint, "-", threshold),
                                             paste0(">", threshold)),
                                  right = T)]
levels(nldas_mle_annual$mle_cat) <- c(levels(nldas_mle_annual$mle_cat),"0")
nldas_mle_annual[is.na(mle_cat), mle_cat := "0"]
nldas_mle_annual[,mle_cat := relevel(mle_cat,ref = "0")]

# keep NLDAS cells that have sufficient trapping data
nldas_polygons <- nldas_polygons %>% 
  arrange(nldasID) %>%
  filter(nldasID %in% nldas_mle_annual$nldasID)

nldas_mle_annual_sf <- right_join(nldas_polygons, nldas_mle_annual, on = "nldasID")

st_bbox(nldas_polygons)

#######################         Training Data        ##########################

# Get training data

# set time bounds for training data
start_year <- 2006
end_year <- 2021

# sort NLDAS cells by month
nldas_ids <- unique(nldas_dt$nldasID)
nldas_dt <- nldas_dt[order(nldasID, year, month)]

# original and new column names for lag
og_cols <- c("temp", "evpsfc")
new_cols <- c("tmp_lag", "evp_lag")

# lag environmental variables by 4 months
# each year will now have data from September to August instead of January to December
env_dt <- nldas_dt[, (new_cols) := shift(.SD, 4), by = nldasID, .SDcols = og_cols]

# calculate standardized environmental values
env_dt_sub <- env_dt[year %in% start_year:end_year]
env_std <- env_dt_sub[,.(year,tmp_std = scale(tmp_lag, scale = T, center = T)[,1],
                         evp_std = scale(evp_lag, scale = T, center = T)[,1]),.(nldasID,month)]

# calculate mean and sd for environmental variables
# used in scaling test data
env_means <- env_dt[year %in% start_year:end_year,.(mean_tmp = mean(tmp_lag,na.rm = T), sd_tmp = sd(tmp_lag,na.rm = T), 
                                                    mean_evp = mean(evp_lag,na.rm = T), sd_evp = sd(evp_lag,na.rm = T)),.(nldasID, month)]

# pivot standardized environmental vars from long to wide
env_wide <- dcast(env_std, nldasID + year ~ month, value.var = c("tmp_std","evp_std"))
setnames(env_wide, c("nldasID","year", paste0("tmp_",c(9:12,1:8)), paste0("evp_",c(9:12,1:8))))

# join wide standardized environmental data to MLE
train_dt <- env_wide[nldas_mle_annual, on = c("nldasID","year"),nomatch = 0]

# generate all combos of env variables for modeling
var.dt <- as.data.table(t(combn(names(train_dt)[3:26], n_preds)))

#######################       Yearly Predictions     ###########################

# Yearly predictions

# get_forecast_for_year
#
# Inputs: 
#   1. year_i -- Year to generate forecasts for
#   2. m.wts -- Parameters and adjusted weights of models in ensemble
#
# Outputs:
#   1. nldas_4p_ens_model_fcst -- A data table containing the following for each NLDAS cell:
#     a. Linear ensemble estimate and variance
#     b. Log estimate and standard error

get_forecast_for_year <- function(year_i, m.wts) {
  
  # find models in ensemble
  ens_models_list <- unique(m.wts$model_num)
  
  # create test dataset standardized by the training means and standard deviations
  test_dt <- env_means[env_dt[year == year_i], on = c("nldasID","month"), nomatch = 0]
  test.std <- test_dt[,.(tmp_std = (tmp_lag - mean_tmp)/sd_tmp, evp_std = (evp_lag - mean_evp)/sd_evp),.(nldasID,month,year)]
  test_dt <- dcast(test.std, nldasID + year ~ month, value.var = c("tmp_std","evp_std"))
  setnames(test_dt, c("nldasID","year", paste0("tmp_",c(9:12,1:8)), paste0("evp_",c(9:12,1:8))))
  
  # get predicted values for models in ensemble
  nldas_4p_all_model_fcst <- rbindlist(lapply(ens_models_list, 
                                              prediction_from_model, 
                                              train_dt = train_dt,
                                              test_dt = test_dt,
                                              var.dt = var.dt,
                                              m.wts = m.wts), idcol = "model_num")
  
  nldas_4p_all_model_fcst[,model_num := ens_models_list[model_num]]
  
  # combine model estimates using weighted sum, calculate estimate variance
  nldas_4p_ens_model_fcst <- unique(nldas_4p_all_model_fcst[,.(est,var,weight,ensemble_est = sum(weight*est)),.(nldasID)][
    ,.(ensemble_est,ensemble_var = (sum(weight*sqrt(var + (ensemble_est - est)^2)))^2),nldasID])
  
  # convert estimate and variance to linear scale
  nldas_4p_ens_model_fcst[,`:=`(mu = log(ensemble_est^2/sqrt(ensemble_var + ensemble_est^2)),
                                sigma = sqrt(log((ensemble_var/ensemble_est^2) + 1)))]
  
  # calculate confidence interval
  nldas_4p_ens_model_fcst[,`:=`(L.95 = exp(mu - (sigma*1.96)),
                                U.95 = exp(mu + (sigma*1.96)))]
  
  nldas_4p_ens_model_fcst[,`:=`(model_num = "Ensemble",
                                year = year_i)]
  
  return(nldas_4p_ens_model_fcst)
  
}

# get Jan-Aug models and weights
model_wts <- get_model_weights(model_ests, month_range = c(11, 7))

# get forecasts for 2022 based on 2006-2021 training data
if (is.character(model_wts)) {
  fcst_2022 <- "No forecast."
} else {
  fcst_2022 <- get_forecast_for_year(year = 2022, m.wts = model_wts)
}

#######################    Yearly Prediction Plots   ##########################

# Observed-Predicted & Agreement Plots

plot_observed_predicted_agreement <- function() {
  
  # get observed MLE for 2022
  obs_mle_2022 <- nldas_mle_annual[year == 2022]
  obs_mle_2022[,model_num := "Observed"]
  
  
  # compare forecasts to observed to determine agreement
  agreement_2022 <- fcst_2022[obs_mle_2022, on = c("nldasID","year")]
  agreement_2022[MLE >= threshold, Agreement := ensemble_est >= threshold]
  agreement_2022[MLE < threshold, Agreement := ensemble_est < threshold]
  
  # prepare agreement values for plotting with observed and predicted
  agreement_2022 <- setDT(melt(agreement_2022, id.vars = c("nldasID"),
                               measure.vars = "Agreement",
                               variable.name = "Prediction", value.name = "MLE", variable.factor = F))
  agreement_2022[,Est.Type := "Agreement"]

  # combine observed and forecast
  obs_fcst <- rbindlist(list(observed = obs_mle_2022[,.(nldasID,MLE,MLE_L95,MLE_U95)] %>% 
                               setnames(c("MLE","MLE_L95","MLE_U95"),c("Est","L.95","U.95")),
                             forecast = fcst_2022[,.(nldasID,ensemble_est,L.95,U.95)] %>% 
                               setnames("ensemble_est","Est")),idcol = "Est.Type")
  
  # create binned MLE variable
  obs_fcst[, MLE := cut(Est, breaks = c(0, midpoint, threshold, 1E30),
                        labels = c(paste0("<", midpoint),
                                   paste0(midpoint, "-", threshold),
                                   paste0(">", threshold)),
                        right = T)]
  levels(obs_fcst$MLE) <- c(levels(obs_fcst$MLE),"0")
  obs_fcst[is.na(MLE), MLE := "0"]
  obs_fcst[,MLE := relevel(MLE,ref = "0")]
  
  # combine observed, predicted, and agreement values
  obs_pred_agreement <- rbindlist(list(obs_fcst, agreement_2022), fill = TRUE)
  
  # change variable types for plotting
  obs_pred_agreement[, MLE := as.character(MLE)]
  levels(obs_pred_agreement$MLE) <- c(levels(obs_fcst$MLE),"TRUE","FALSE")
  obs_pred_agreement[, Est.Type := case_when(Est.Type == "observed" ~ "Observed",
                                             Est.Type == "forecast" ~ "Predicted",
                                             Est.Type == "Agreement" ~ "Agreement")]
  obs_pred_agreement[,Est.Type := factor(Est.Type, levels = c("Observed","Predicted","Agreement"))]
  
  # combine infection data with NLDAS cell polygons
  obs_pred_agreement_sf <- left_join(nldas_polygons, obs_pred_agreement, by = "nldasID")
  obs_pred_agreement_sf$MLE <- factor(obs_pred_agreement_sf$MLE, levels = levels(obs_pred_agreement$MLE))
  
  # colors for plotting
  pred_cols <- c(RColorBrewer::brewer.pal(name = "YlOrRd", n = 4), "#06ABEB", "#DC298D")
  names(pred_cols) <- levels(obs_pred_agreement$MLE)
  
  # plot observed, predicted, and agreement values together
  pred_map_2022 <- ggmap(terrain) + 
    geom_sf(data = nldas_polygons, fill = NA, inherit.aes = F) + 
    geom_sf(data = obs_pred_agreement_sf %>% filter(!is.na(MLE)), aes(fill = MLE), inherit.aes = F) + 
    scale_fill_manual(values = pred_cols) + 
    geom_sf_text(data = obs_pred_agreement_sf %>% filter(!is.na(MLE) & Est.Type != "Agreement"), 
                 aes(label = round(Est,1)), size = 6, inherit.aes = F) + 
    facet_grid(~ Est.Type, switch = "y") + 
    labs(y = "Latitude",
         x = "Longitude",
         fill = expression(I[M]), 
         title = "2022 Observed and Ensemble-Predicted Infection Rates with Agreement") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.subtitle = element_blank(),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))
  
  ggsave(paste0(figure_folder, "obs_pred_agreement.png"), pred_map_2022, width = 12, height = 5)
  
}

plot_observed_predicted_agreement()

#######################      Monthly Predictions     ##########################

# 2022 Monthly Predictions

get_monthly_predictions <- function(pred_month, m.wts) {
  
  # find models in ensemble
  ens_models_list <- unique(m.wts$model_num)
  
  # prior_dt is environmental data for months prior to forecast date
  # environmental data after the forecast date is set as the mean from previous years
  # this standardizes to 0 in the next step
  prior_dt <- env_dt[year == 2022 & month %in% 1:(pred_month + 4), .(nldasID, month, tmp_lag, evp_lag)]
  future_dt <- env_means[month %in% (pred_month + 5):12, .(nldasID, month, tmp_lag = mean_tmp, evp_lag = mean_evp)]
  raw_test_dt <- setDT(rbind(prior_dt, future_dt))
  
  # create test data
  test_dt <- env_means[raw_test_dt, on = c("nldasID", "month")]
  
  # standardize test data using prior environmental data means and standard deviations
  test.std <- test_dt[,.(tmp_std = (tmp_lag - mean_tmp)/sd_tmp, evp_std = (evp_lag - mean_evp)/sd_evp),.(nldasID,month)]
  test_dt <- dcast(test.std, nldasID ~ month, value.var = c("tmp_std","evp_std"))
  setnames(test_dt, c("nldasID", paste0("tmp_",c(9:12,1:8)), paste0("evp_",c(9:12,1:8)))) 
  
  # generate forecast for each model in ensemble
  nldas_4p_all_model_fcst <- rbindlist(lapply(ens_models_list, 
                                              prediction_from_model, 
                                              train_dt = train_dt,
                                              test_dt = test_dt, 
                                              var.dt = var.dt,
                                              m.wts = m.wts), idcol = "model_num")
  nldas_4p_all_model_fcst[,model_num := ens_models_list[model_num]]
  
  # combine model estimates using a weighted sum, calculate variance of estimate
  nldas_4p_ens_model_fcst <- unique(nldas_4p_all_model_fcst[,.(est,var,weight,ensemble_est = sum(weight*est)),.(nldasID)][
    ,.(ensemble_est,ensemble_var = (sum(weight*sqrt(var + (ensemble_est - est)^2)))^2),nldasID])
  
  # convert from log-scale to linear scale
  nldas_4p_ens_model_fcst[,`:=`(mu = log(ensemble_est^2/sqrt(ensemble_var + ensemble_est^2)),
                                sigma = sqrt(log((ensemble_var/ensemble_est^2) + 1)))]
  
  # calculate confidence interval
  nldas_4p_ens_model_fcst[,`:=`(L.95 = exp(mu - (sigma*1.96)),
                                U.95 = exp(mu + (sigma*1.96)))]
  
  nldas_4p_ens_model_fcst[,`:=`(model_num = "Ensemble",
                                month = pred_month)]
  
  return(nldas_4p_ens_model_fcst)
  
}

# get Nov-Aug models and weights
model_wts <- get_model_weights(model_ests, month_range = c(11, 7))

if (is.character(model_wts)) {
  forecasts_mar_jul_2022 <- "No forecast."
} else {
  forecasts_mar_jul_2022 <- lapply(3:7, get_monthly_predictions, m.wts = model_wts) %>% bind_rows()
}

#######################    Monthly Prediction Plots  #########################

# 2022 Monthly Observed-Predicted & Agreement Plots

plot_observed_predicted_agreement_monthly <- function() {
  
  if(is.character(forecasts_mar_jul_2022)) {
    return("The ensemble was not created. No plots could be generated.")
  }
  
  # get observed MLE for 2022
  obs_mle_2022 <- nldas_mle_annual[year == 2022]
  obs_mle_2022[,model_num := "Observed"]
  
  obs_mle_sf <- left_join(nldas_polygons, obs_mle_2022)

  # compare forecasts to observed to determine agreement
  agreement_2022 <- forecasts_mar_jul_2022[obs_mle_2022, on = c("nldasID")]
  agreement_2022[MLE >= 1,Agreement := ensemble_est >= 1]
  agreement_2022[MLE < 1,Agreement := ensemble_est < 1]
  
  # prepare agreement values for plotting with forecast
  agreement_2022 <- setDT(melt(agreement_2022, id.vars = c("nldasID", "month"),
                               measure.vars = "Agreement",
                               variable.name = "Prediction", value.name = "MLE", variable.factor = F))
  agreement_2022[,Est.Type := "Agreement"]

  # prepare forecast values for plotting
  obs_fcst <- forecasts_mar_jul_2022[,.(nldasID,ensemble_est,L.95,U.95, month)] %>%
    setnames("ensemble_est","Est") %>%
    mutate(Est.Type = "Predicted") %>%
    setDT()

  # create binned MLE variable
  obs_fcst[, MLE := cut(Est, breaks = c(0, midpoint, threshold, 1E30),
                        labels = c(paste0("<", midpoint),
                                   paste0(midpoint, "-", threshold),
                                   paste0(">", threshold)),
                        right = T)]
  levels(obs_fcst$MLE) <- c(levels(obs_fcst$MLE),"0")
  obs_fcst[is.na(MLE), MLE := "0"]
  obs_fcst[,MLE := relevel(MLE,ref = "0")]

  # combine observed, predicted, and agreement values
  obs_pred_agreement <- rbindlist(list(obs_fcst, agreement_2022), fill = TRUE)
  
  # change variable types for plotting
  obs_pred_agreement[, MLE := as.character(MLE)]
  levels(obs_pred_agreement$MLE) <- c(levels(obs_fcst$MLE),"TRUE","FALSE")
  obs_pred_agreement[,Est.Type := factor(Est.Type, levels = c("Predicted","Agreement"))]

  # combine infection data with NLDAS cell polygons
  obs_pred_agreement_sf <- left_join(nldas_polygons, obs_pred_agreement, by = "nldasID")
  obs_pred_agreement_sf$MLE <- factor(obs_pred_agreement_sf$MLE, levels = levels(obs_pred_agreement$MLE))

  # set plot colors
  pred_cols <- c(RColorBrewer::brewer.pal(name = "YlOrRd", n = 4), "#06ABEB", "#DC298D")
  names(pred_cols) <- levels(obs_pred_agreement$MLE)
  
  # plot observed infection rates for 2022
  observed <- ggmap(terrain) +  
    geom_sf(data = region_shp, fill = NA, inherit.aes = F, colour = "black", size = 2) + 
    geom_sf(data = nldas_polygons, fill = NA, inherit.aes = F) + 
    geom_sf(data = obs_mle_sf, aes(fill = mle_cat), inherit.aes = F) +
    geom_sf_text(data = obs_mle_sf, aes(label = round(MLE,1)), size = 7, inherit.aes = F) + 
    scale_fill_manual(values = pred_cols) +
    labs(x = "Longitude", 
         y = "Latitude",
         title = NULL,
         fill = expression(I[M])) + 
    facet_grid(~model_num) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin=grid::unit(c(0, 4, 0, 2), "cm"),
          axis.ticks.y = element_blank(),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    guides(fill = guide_legend(override.aes = list(size = 0.5)))
  
  ggsave(paste0(figure_folder, "observed_2022.png"), observed, width = 8, height = 8)
  
  # set labels for monthly predictions plot
  month_names <- paste(month.name[c(4:8)], 4)
  names(month_names) <- as.character(c(3:7))
  
  # plot monthly predictions and agreement
  pred_map_2022 <- ggmap(terrain) + 
    geom_sf(data = nldas_polygons, fill = NA, inherit.aes = F) + 
    geom_sf(data = obs_pred_agreement_sf %>% filter(!is.na(MLE)), aes(fill = MLE), inherit.aes = F) + 
    scale_fill_manual(values = pred_cols) + 
    geom_sf_text(data = obs_pred_agreement_sf %>% filter(!is.na(MLE) & Est.Type != "Agreement"), 
                 aes(label = round(Est,1)), size = 3, inherit.aes = F) + 
    facet_grid(Est.Type ~ month, 
               labeller = labeller(month = month_names)) + 
    labs(y = "Latitude",
         x = "Longitude",
         fill = expression(I[M])) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.subtitle = element_blank(),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))
  
  ggsave(paste0(figure_folder, "pred_agreement_monthly.png"), pred_map_2022, width = 12, height = 6)
  
}

plot_observed_predicted_agreement_monthly()

#######################       Scatterpie Plots       ###########################

# Scatterpie plotting

# identify environmental variables in model (temperature and evapotranspiration)
vars <- c("tmp", "evp")
var_labels <- c("TEMP", "ET")

# get_rel_imp
#
# Inputs:
#   1. m.wts -- The weights of the model subset (Jan-Aug or Nov-Aug)
#   2. model_ests -- All model parameters
#   3. nov -- Boolean representing if subset is Nov-Aug (TRUE = Nov-Aug)
#
# Outputs:
#   1. Relative importance of each variable in the ensemble by month

get_rel_imp <- function(model_wts, start_month) {
  
  # remove intercept, extract sign from each variable's estimate
  month_wts <- model_wts[varnames != "(Intercept)",.(dir = sign(Estimate),adj.wts),.(varnames,model_num)]
  
  # weighted proportion of each variable in ensemble grouped by direction
  relative_importance <- unique(month_wts[,.(relative_imp = sum(adj.wts)), .(dir, varnames)])
  relative_importance[,c("variable","month") := tstrsplit(varnames,"_")]
  relative_importance[,month := as.numeric(month)]
  
  # split env variables
  
  for (i in 1:length(vars)) {
    
    relative_importance[variable == vars[i], var_position := i]
    
  }
  
  # pivot relative_importance into negative and positive portions
  ri.4p <- dcast(relative_importance,variable + month + var_position ~ dir,value.var = "relative_imp", fill = 0)
  setnames(ri.4p,c('-1','1'),c("neg_perc","pos_perc"))
  
  # placeholder to fill rest of pie plot
  
  ri.4p[,no_wt := 1 - (neg_perc + pos_perc)]
  
  # shift month labels for November and December plotting
  if(start_month >= 9) {
    ri.4p[, shifted_month := ifelse(month > 8, month - start_month + 1, month + 13 - start_month)]
  } else {
    ri.4p[, shifted_month := (month - start_month + 1)]
  }
  
  # return relative importance data table
  return(ri.4p)
  
}

# plot_scatterpie
#
# Inputs:
#   1. model_ests -- All model parameters
#   2. subtitle -- Subtitle for the scatterpie plot, use default as guide
#
# No direct outputs, saves scatterpie plot for November to July

plot_scatterpie <- function(model_ests,
                            subtitle = "NLDAS scale 4-predictor ensemble model, 2006 to 2022") {
  
  # define positive, negative, and filler colors, respectively
  pos_cols <- c("#06ABEB","#DC298D", "white")
  
  # get ensemble weights for variables between November and July, inclusive
  nov_jul_wts <- get_model_weights(model_ests, c(11, 7))
  
  if (!is.character(nov_jul_wts)) {
    
    write.csv(nov_jul_wts, paste0(figure_folder, "nov_jul_models_", start_year, "_", end_year, ".csv"), row.names = FALSE)
    nov_jul_num_models <- nrow(nov_jul_wts) / 5
    
    # get relative importance of each variable in subsetted ensemble
    nov_jul_ri <- get_rel_imp(nov_jul_wts, start_month = 11)
    
    # Nov-Jul scatterpie
    
    nov_jul_nldas_4p_subset <- ggplot() +
      geom_scatterpie(data = nov_jul_ri,
                      aes(x = shifted_month, y = var_position),
                      pie_scale = 3,
                      cols = c("pos_perc","neg_perc","no_wt")) +
      coord_fixed() +
      scale_x_discrete(limits = c(month.abb[c(11:12, 1:7)]), labels = c(month.abb[c(11:12, 1:7)])) +
      scale_y_continuous(breaks = 1:length(vars), labels = var_labels) +
      scale_fill_manual(values = pos_cols, 
                        name = expression("Effect on" ~ I[M]),
                        labels = c("Increase", "Decrease", ""),
                        breaks = c("pos_perc", "neg_perc", "no_wt")) +
      labs(x = "Month",title = "Environmental Variable Importance by Month, Coachella Valley",
           subtitle = subtitle,
           caption = paste0("\n\nNote: ", nov_jul_num_models, " models with environmental variables in Nov-Jul included in ensemble model")) +
      ylab(NULL) +
      theme(plot.title = element_text(size = 18, face = "bold"),
            plot.subtitle = element_text(size = 16),
            plot.caption = element_text(size = 12, hjust = 0))
    
    nov_jul_nldas_4p_subset
    
    # save plots to specified destinations
    
    ggsave(nov_jul_nldas_4p_subset,
           filename = paste0(figure_folder, "CV_NLDAS_4p_", start_year, "_", end_year, "_Nov_Jul.png"),
           width = 15, height = 8)
    
  }
  
}

plot_scatterpie(model_ests,
                subtitle = paste("NLDAS scale 4-predictor ensemble model, ", 
                                 start_year, " to ", end_year, "\n", MLE_def_str, 
                                 " \n", "\n", sep = ""))

#######################    Annual Infection Rates    ###########################

# Annual infection rates plotting

annual_inf_rates <- ggmap(terrain) + 
  geom_sf(data = region_shp, inherit.aes = F, fill = NA) + 
  geom_sf(data = nldas_polygons, inherit.aes = F, fill = NA, colour = "light grey") + 
  geom_sf(data = nldas_polygons, inherit.aes = F, fill = NA, colour = "black") + 
  geom_sf(data = nldas_mle_annual_sf,inherit.aes = F, 
          aes(fill = mle_cat)) + 
  geom_sf_text(data = nldas_mle_annual_sf,
               aes(label = round(MLE,1)), size = 3, inherit.aes = F) +
  scale_fill_brewer(palette = "YlOrRd", type = "seq") + 
  # change axis labels and title
  labs(x = "Longitude", y ="Latitude", 
       title = "Annual WNV Mosquito Infection Rate, Coachella Valley",
       fill = expression(I[M])) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10),
        legend.position = c(0.92, 0.1),
        legend.justification = c(0.9, 0.1),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) + 
  facet_wrap(~year, nrow = 3, ncol = 6)

ggsave(annual_inf_rates, 
       filename = paste0(figure_folder, "NLDAS_Annual_Inf_Rates.png"),
       height = 10, width = 16, units = "in")

#######################    NLDAS ID + Frequencies    ###########################

# NLDAS ID Plotting

nldas_id_counts <- nldas_mle_annual_sf %>%
  group_by(nldasID) %>%
  summarize(count = n())

nldas_id_map <- ggmap(terrain) +  
  geom_sf(data = region_shp, fill = NA, inherit.aes = F, colour = "black", size = 2) + 
  geom_sf(data = nldas_polygons, fill = "white", inherit.aes = F) + 
  geom_sf_text(data = nldas_id_counts, aes(label = paste(nldasID, "\n", count, sep = "")), size = 3.5, inherit.aes = F) + 
  labs(x = NULL, y = NULL, title = "NLDAS IDs") + 
  theme(axis.text = element_blank(),
        plot.margin=grid::unit(c(1,1,1,0.5), "cm"),
        axis.ticks.y = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 0.5)))

ggsave(paste0(figure_folder, "nldas_id_map.png"), nldas_id_map, width = 8, height = 8)


#######################        MLE Histogram         ###########################

# histogram to show distribution of infection rates in training data
mle_hist <- nldas_mle_annual %>%
  filter(year %in% start_year:end_year) %>%
  ggplot(aes(x = MLE)) +
  geom_histogram() + 
  geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
  annotate("text", x = threshold + 0.2, y = 60, label = "75th Percentile", angle = 90) + 
  scale_x_continuous(name = "Mean Annual Mosquito Infection Rate Per 1000 Tested",
                     breaks = seq(0, 30, 5)) + 
  labs(y = "Frequency")

ggsave(paste0(figure_folder, "mle_histogram.png"), mle_hist,
       width = 10, height = 6)

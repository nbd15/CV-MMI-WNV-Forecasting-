source("glmm_functions.R")

library(writexl)

#######################################################################################
set.seed(153)

# import NLDAS and MLE data
nldas_dt <- fread("data/cv_monthly_nldas_2005_2022.csv")
nldas_mle_annual <- fread("data/cv_n500_mle_2006_2022_over_11_obs.csv")

# convert NLDAS IDs to factor and find all cells with environmental data
nldas_dt[,nldasID := factor(nldasID)]
nldas_mle_annual[,nldasID := factor(nldasID)]
nldas_ids <- unique(nldas_dt$nldasID)

nldas_dt <- nldas_dt[order(nldasID, year, month)]

# set number of cores for parallel processing
ncores = 32

start_year <- 2006
end_year <- 2021

years <- start_year:end_year

# 2006 - 2021 Ensemble model

# original and new column names for lagging
og_cols <- c("temp", "evpsfc")
new_cols <- c("tmp_lag", "evp_lag")

# shift environmental data back by 4 months
# each year is now represented as September to August
env_dt <- nldas_dt[, (new_cols) := data.table::shift(.SD, 4), by = nldasID, .SDcols = og_cols]
env_dt <- env_dt[year %in% years]

# standardize training data
env_std <- env_dt[,.(year,tmp_std = scale(tmp_lag, scale = T, center = T)[,1],
                     evp_std = scale(evp_lag, scale = T, center = T)[,1]),.(nldasID,month)]

# pivot from long to wide
env_wide <- dcast(env_std, nldasID + year ~ month, value.var = c("tmp_std","evp_std"))
setnames(env_wide, c("nldasID","year",
                     paste0("tmp_",c(9:12,1:8)),
                     paste0("evp_",c(9:12,1:8))))

# join MLE data with environmental data
train_dt <- env_wide[nldas_mle_annual, on = c("nldasID","year"),nomatch = 0]

## 4 predictors

n_preds <- 4

# find all unique combinations of environmental variables
var.dt <- as.data.table(t(combn(names(train_dt)[3:26], n_preds)))

# 2) Run hierarchical model for each row of the dataset above
# mclapply splits list of models amongst 32 cores to speed up processing
dt_glmm_nbinom <- mclapply(1:nrow(var.dt),function(row_num){
  
  tryCatch(admb_fun(i = row_num,model_name = "nbinom", 
                    # dt is the input dataset with the env. vars and inf. rates for each grid cell
                    # dt_var is the dataset with all combinations of the environmental variables
                    # i is the row number of the dt_var dataset
                    # spat_var is the name of the spatial variable to group by
                    # model_name is the name of the family to model the residuals
                    dt = train_dt,spat_var = "nldasID",dt_var = var.dt),
           error=function(e){cat("ERROR :",conditionMessage(e), "\n")})},mc.cores = ncores)

# Identify which set of variables failed to converge and set these to NA
failed <- which(sapply(dt_glmm_nbinom, is.null))

dt_glmm_nbinom[failed] <- lapply(1:length(failed), function(void){
  data.table(Estimate = NA_real_, `Std. Error` = NA_real_, 
             `z value` = NA_real_, pval = NA_real_, 
             varnames = NA_real_, aicc = NA_real_)})

# 3) Get a dataset of coefficients, variables, p-values and AICc score for each model that converged
dt_glmm_nbinom_df <- rbindlist(dt_glmm_nbinom,idcol = T)

setnames(dt_glmm_nbinom_df,".id","model_num")

# Calculate Akaike weights for each model in the dataset from Step 3
dt_sigModels <- glmm_wts(dt_glmm_nbinom_df)

output = list("Model Parameters" = dt_glmm_nbinom_df,
              "Model Weights" = dt_sigModels)

write_xlsx(output, paste0("model_outputs/", start_year, "_", end_year, "/cv_", start_year, "_", end_year, "_all_years_over_11_obs.xlsx"))

suppressPackageStartupMessages({
  library(lubridate)
  library(data.table)
  library(stringr)
  library(tidyverse)
  library(glmmADMB)
  library(glmulti)
  library(MASS)
  library(parallel)
})

#####################################################################################

# Functions to use

admb_fun <- function(i, model_name, dt,spat_var,dt_var,...){
  
  # Arguments: 
  # dt is the input dataset with the env. vars and inf. rates for each grid cell
  # dt_var is the dataset with all combinations of the environmental variables
  # i is the row number of the dt_var dataset
  # spat_var is the name of the spatial variable to group by
  # model_name is the name of the family to model the residuals
  
  # This model selects the combination of variables present in row i of the 
  # dataset dt_var, subsets the original dataset to these variables along with the 
  # dependent variable and grouping factor, and runs a 
  # hierarchical model for count data, with a random intercept for the 
  # grouping variable and a residual error family specified according to the 
  # model_name argument. 
  # The function outputs a dataframe of the coefficents of the regression model,
  # the SEs of the coefficients, p-values from the Wald tests of the coefficient and
  # the AICc of the model.
  
  message("Working with row ",i)
  
  # Subset input dataset
  test.subset <- dt[,names(dt)[match(c("MLE",spat_var,dt_var[i,]), 
                                     names(dt))],with = F]
  
  # Build the random effects
  rand_term <- paste0("(",1,"|",spat_var,")")
  
  # build the regression formula
  f = as.formula(paste0("MLE~",paste(paste(dt_var[i,],collapse = "+"),
                                     rand_term, sep = "+")))
  
  if(model_name == "poisson"){
    admb.model <- glmmadmb(f, data = test.subset, family = "poisson",zeroInflation = T)
  }else if(model_name == "nbinom"){
    admb.model <- glmmadmb(f, data = test.subset, family = "nbinom",zeroInflation = T)
  } else if (model_name == "nbinom1"){
    admb.model <- glmmadmb(f, data = test.subset, family = "nbinom1",zeroInflation = T)
  }
  
  out <- summary(admb.model)
  t <- as.data.table(out$coef)
  t[,varnames := row.names(out$coef)]
  t[,aicc := aicc(admb.model)]
  setnames(t, old = "Pr(>|z|)", new = "pval")
  return(t)
}

glmm_wts <- function(dt,p.val = 0.05,...){
  
  # Function to calculate Akaike weights
  # Given a dataset 'dt' consisting of model variables, coefficients, p-values and AICc 
  # for each model, this function identifies the relative weights that each model carries
  
  # Identify and report the number of terms that are significant in each model, and the largest such number.
  max_sig_terms <- dt[!is.na(pval),.(n_terms = .N,
                                     sig_terms = sum(pval<= p.val),
                                     aicc),.(model_num)][,max(sig_terms)]
  
  
  if (max_sig_terms != dt[!is.na(pval), .(n_terms = .N), .(model_num)][,max(n_terms)]) {
    
    message("The ensemble could not be created. No models were fully statistically significant.")
    return("No significant models")
    
  }
  
  # Select all models where the number of significant terms is equal to the maximum number of sig. terms,
  # and order these in ascending order of AICc score. 
  sigModels <- unique(dt[!is.na(pval),.(n_terms = .N,
                                        sig_terms = sum(pval<= p.val),
                                        aicc),.(model_num)][sig_terms == max_sig_terms])[order(aicc)]
  
  # Calculate difference between AICc score of each model and AICc score of best model.
  sigModels[,delta := aicc - min(aicc)]
  
  # Exponentiate the difference
  sigModels[,exp_delta := exp((-1/2)*delta)]
  
  # Normalize exponentiated difference over all models - the Akaike weight of each model 
  sigModels[,weights := exp_delta/sum(exp_delta)]
  
  # Get the cumulative Akaike weight over all the models
  sigModels[,cum.wts := cumsum(weights)]
  return(sigModels)
}

get_model_weights <- function(model_ests, month_range = c(11, 7), vars = c("tmp", "evp"), n_preds = 4) {
  
  # get_model_weights
  #
  # Inputs:
  #   1. model_ests -- All model parameters
  #   2. month_range -- Range of months to include in ensemble
  #     a. Integer vector with default c(11, 7): November to July
  #
  # Outputs:
  #   1. model_wts -- A data table with the parameters and adjusted weights for each model in the ensemble
  
  start_month <- month_range[1]
  end_month <- month_range[2]
  inc_months <- c()
  
  # find months to be included in the ensemble
  if (start_month >= 9) {
    
    inc_months <- c(start_month:12, 1:end_month)
    
  } else inc_months <- start_month:end_month
  
  # find variables within specified month range
  vars <- as.vector(unlist(sapply(paste0(vars, "_"),
                                  function(x) paste0(x,inc_months), simplify = F)))
  
  # find models that only contain variables within month range
  models_month_subset <- model_ests[,sum(varnames %in% vars), model_num][V1 == n_preds, unique(model_num)]
  message("Calculating ensemble fixed effects and ensemble variance")
  
  model_ests <- model_ests[model_num %in% models_month_subset]
  
  # calculate ensemble weights for model subset
  # glmm_wts found in glmm_functions.R
  m.wts <- glmm_wts(model_ests)
  
  if(is.character(m.wts)) {
    return(m.wts)
  }
  
  # find the models that will be in the final ensemble
  # the first model with >95% cumulative weight will be the last model included
  m.wts[, row.num := .I]
  ens_model_rn <- m.wts[cum.wts > 0.95, row.num[1]]
  ens_models_list = m.wts[1:ens_model_rn, model_num]
  m.wts <- m.wts[model_num %in% ens_models_list]
  
  # adjust weights based on total weights of model subset
  m.wts[, adj.wts := weights/sum(weights)]
  model_wts <- m.wts[, .(model_num, adj.wts, weights)][model_ests, on = "model_num", nomatch = 0]
  
  # return data table with each model's output and its adjusted weight in the ensemble
  return(model_wts)
  
}

prediction_from_model <- function(i.model, train_dt, test_dt, var.dt, m.wts){
  
  # prediction_from_model
  #
  # Inputs:
  #   1. i.model -- Model number to generate prediction for
  #   2. test_dt -- Test environmental data to get prediction for
  #   3. m.wts -- Parameters and adjusted weights of models in ensemble
  #
  # Outputs:
  #   1. final.lin.est -- A data table containing the following for each NLDAS cell:
  #     a. Fixed and random effects
  #     b. Estimates in linear and log scale
  #     c. Variance and confidence intervals in linear scale
  #     d. Weight associated with the model used for predictions
  
  message(i.model)
  
  # Get the re-normalised weight for the model selected
  wt <- unique(m.wts[model_num == i.model, adj.wts])
  
  # Get the data for the variables in the model
  test.subset <- train_dt[,names(train_dt)[match(c("MLE","nldasID",var.dt[i.model,]), 
                                                 names(train_dt))],with = F]
  
  f = as.formula(paste0("MLE~",paste(paste(var.dt[i.model,],collapse = "+"),
                                     "(1|nldasID)", sep = "+")))
  # Run the model
  admb.model <- glmmadmb(f, data = test.subset, family = "nbinom",zeroInflation = T)
  
  # Ensure the test and train datasets have the same nldasIDs, otherwise prediction throws an error
  
  nids <- intersect(train_dt$nldasID, test_dt$nldasID)
  
  test.dt.subset <- test_dt[nldasID %in% nids]
  
  # Get the fixed effects and their standard errors from the model for the IDs in nids above, on the log scale
  FEs <- predict(admb.model, type = "link",
                 newdata = test.dt.subset, se.fit = T)
  REs <- ranef(admb.model)$nldasID # Get the random effects for each nldasID in the test dataset subset
  
  # REs in the step above are a matrix, with the row names being the nldasIDs. Attache these row names to enable easier
  # identification
  REs <- REs[rownames(REs) %in% nids,]
  
  # Get the prediction from the model on the log scale as the sum of the fixed effect and random effect for each 
  # NLDAS grid cell 
  final.lin.est <- as.data.table(cbind(fixed = FEs$fit,random = REs,
                                       log.est = REs + FEs$fit, # log scale estimate for each NLDAS grid
                                       log.se = FEs$se.fit # log scale SE for each NLDAS grid
  ),keep.rownames = "nldasID")    
  
  # I am calculating the confidence intervals on the natural scale in this step first, 
  # even though the estimates and standard errors are on the log scale. 
  # I calculate the estimates on the natural scale further below within the function.
  final.lin.est[,`:=`(L.95 = exp(log.est - (log.se*1.96)),
                      U.95 = exp(log.est + (log.se*1.96))),nldasID]
  
  # Convert the prediction on the log scale to the natural scale, with the appropriate adjustment from the log-normal 
  # distribution to the normal distribution
  final.lin.est[,est := exp(log.est + (log.se^2)/2)]
  final.lin.est[,var := exp((2*log.est) + log.se^2)*(exp(log.se^2) - 1)]
  final.lin.est[,weight := wt]
  
  return(final.lin.est)
}

nldas_prediction_function <- function(train.dt,test.dt,n_predictors,
                                      varfrom,varto,ncores,month_range,...){
  
  # NLDAS Prediction Function
  # This function outputs predicted infection rate from a hierarchical neg-binom/poisson model.
  # The temporal scale of the prediction depends on the inputs to the model;
  # the spatial scale is the NLDAS grid.
  
  # Description of arguments: 
  # train.dt = all years of data except year in target
  # test.dt = data for year in target
  # n_predictors = number of env. variables in the model
  # varfrom, varto = identify the index number of env. variables in the training dataset, which are used to 
  # generate combinations
  
  # 1) Generate dataset with all possible combinations of environmental variables.
  var.dt <- as.data.table(t(combn(names(train.dt)[varfrom:varto],n_predictors)))
  
  # 2) Run hierarchical model for each row of the dataset above
  dt_glmm_nbinom <- mclapply(1:nrow(var.dt),function(row_num){
    
    tryCatch(admb_fun(i = row_num,model_name = "nbinom", 
                      # dt is the input dataset with the env. vars and inf. rates for each grid cell
                      # dt_var is the dataset with all combinations of the environmental variables
                      # i is the row number of the dt_var dataset
                      # spat_var is the name of the spatial variable to group by
                      # model_name is the name of the family to model the residuals
                      dt = train.dt,spat_var = "nldasID",dt_var = var.dt),
             error=function(e){cat("ERROR :",conditionMessage(e), "\n")})},mc.cores = 32) #
  
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
  #dt_sigModels <- glmm_wts(dt_glmm_nbinom_df)
  
  start_month <- month_range[1]
  end_month <- month_range[2]
  
  # find all months that should be included in ensemble
  inc_months <- c()
  
  if (start_month >= 9 & end_month < 9) {
    
    inc_months <- append(inc_months, c(start_month:12, 1:end_month))
    
  } else inc_months <- start_month:end_month
  
  vars <- as.vector(unlist(sapply(c("tmp_","evp_"),
                                  function(x)paste0(x,inc_months), simplify = F)))
  
  models_month_subset <- dt_glmm_nbinom_df[,sum(varnames %in% vars), model_num][V1 == n_predictors, unique(model_num)]
  message("Calculating ensemble fixed effects and ensemble variance")
  
  model_ests <- dt_glmm_nbinom_df[model_num %in% models_month_subset]
  
  dt_sigModels <- glmm_wts(model_ests)
  
  # Get the model numbers and the associated weights that would go into the ensemble
  dt_sigModels[,row.num := .I] # assign row numbers to dataset of weights
  ens_models <- dt_sigModels [cum.wts > 0.95,model_num[1]] # Identify first model number whose cumulative weight is > 0.95
  ens_model_rn <- dt_sigModels[model_num == ens_models, row.num] # get row numbers of model whose cumulative weigh > 0.95
  ens_models_list = dt_sigModels[1:ens_model_rn,model_num] # model numbers whose cumulative weight adds up to 0.95
  # subset dataset of model weights to models whose weights add up to 0.95
  dt_wts.Subset = dt_sigModels[model_num %in% ens_models_list] 
  dt_wts.Subset[,adj.cum.wts := weights/sum(weights)] # re-normalize the cumulative weights
  
  model_subsets <- dt_wts.Subset[,model_num] # get the model numbers in the ensemble model
  best_model <- dt_wts.Subset[,model_num[which.max(weights)]]
  
  # For each model identified in the ensemble, run the hierarchical again with the subset of variables in it,
  # and get predictions
  all_models_out <- rbindlist(lapply(model_subsets, 
                                     prediction_from_model, 
                                     train_dt = train.dt,
                                     test_dt = test.dt,
                                     var.dt = var.dt,
                                     m.wts = dt_wts.Subset), idcol = "model_num")
  
  all_models_out[,model_num := model_subsets[model_num]]
  
  # Once the prediction from each model is available, weight the prediction and 
  # get the ensembled prediction and the SE of the ensembled prediction
  final.ens.est <- unique(all_models_out[,.(est,var,weight,ensemble_est = sum(weight*est)),.(nldasID)][
    ,.(ensemble_est,ensemble_var = (sum(weight*sqrt(var + (ensemble_est - est)^2)))^2),nldasID
  ])
  
  final.ens.est[,`:=`(mu = log(ensemble_est^2/sqrt(ensemble_var + ensemble_est^2)),
                      sigma = sqrt(log((ensemble_var/ensemble_est^2) + 1)))]
  
  final.ens.est[,`:=`(L.95 = exp(mu - (sigma*1.96)),U.95 = exp(mu + (sigma*1.96)))]
  
  # Get prediction for the best model, which is the prediction with the largest weight
  best.est <- all_models_out[model_num == best_model]
  
  return(list(final.est = final.ens.est, best.est = best.est,
              models = dt_glmm_nbinom_df,
              model_weights = dt_sigModels))
}
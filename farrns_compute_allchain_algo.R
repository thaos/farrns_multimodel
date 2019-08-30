source("check_H0_algo.R")
source("CramerVonMisesTwoSamples.R")
compute_allchain <- function(lmodels, tas_cmip5, kernel, bandwidth){
  
  tas_hadcrut <- tas_cmip5[tas_cmip5$model == "HadCRUT", ]
  tas_hadcrut_counterfactual <- tas_hadcrut[tas_hadcrut$year <= 1900,]
  
  lp12 <- lapply(lmodels, function(model, tas_cmip5, kernel, bandwidth){
    print("---------------------------")
    print(model)
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" | tas_model$experiment == "rcp85", ]
    print(head(tas_model_factual))
    if(model == "HadCRUT"){
      tas_model_counterfactual <- tas_hadcrut_counterfactual
    } else {
      tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
    }
    print(head(tas_model_counterfactual))
    # tas_model_counterfactual <- tas_model[tas_model$experiment == "historical" & tas_model$year <= 1900, ]
    tpred <- min(tas_model$year):max(tas_model$year)
    p12 <- estim_p12_ns(x = tas_model_counterfactual$tas, t = tas_model_factual$year, z = tas_model_factual$tas, tpred = tpred, kernel = kernel, h = bandwidth)
  }, tas_cmip5 = tas_cmip5, kernel = kernel, bandwidth = bandwidth)
  
  lUhat <- lapply(lmodels, function(model, tas_cmip5, tas_hadcrut, kernel, bandwidth){
    # browser()
    print("---------------------------")
    print(model)
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical", ]
    if(model == "HadCRUT"){
      tas_model_counterfactual <- tas_hadcrut_counterfactual
    } else {
      tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
    }
    yend <- min(c(max(tas_hadcrut$year), max(tas_model_factual$year), max(tas_model_counterfactual$year)))
    ystart <- max(c(min(tas_hadcrut$year), min(tas_model_factual$year), min(tas_model_counterfactual$year)))
    year_incommon <- intersect(tas_hadcrut$year, tas_model_factual$year) %>%
      intersect(tas_model_counterfactual$year) 
    tas_hadcrut <- tas_hadcrut[tas_hadcrut$year %in% year_incommon, ]
    tas_model_counterfactual <- tas_model_counterfactual[tas_model_counterfactual$year%in% year_incommon, ]
    tas_model_factual <- tas_model_factual[tas_model_factual$year %in% year_incommon, ]
    # if(model == "HadGEM2-ES") browser()
    Uhat <- compute_Uhat(
      t = tas_hadcrut$year,
      tm = tas_model_counterfactual$year,
      Xm = tas_model_counterfactual$tas,
      Z = tas_hadcrut$tas,
      Zm = tas_model_factual$tas,
      kernel = kernel, bandwidth = bandwidth
    )
    return(Uhat)
  }, tas_cmip5 = tas_cmip5, tas_hadcrut = tas_hadcrut, kernel = kernel, bandwidth = bandwidth)
  
  CvM_df <- CvM_checkH0(
    X = tas_hadcrut_counterfactual$tas,
    lUhat = lUhat,
    lmodels = lmodels
  )
  CvM_df_filtered <- keep_onemodel_perinstitute(CvM_df, unique(tas_cmip5[, c("institute", "model")]))
  CvM_df_filtered <- subset(CvM_df_filtered, model != "HadCRUT")
  lmodels_tokeep <- CvM_df_filtered$model
  lp12_tokeep <- lp12[lmodels %in% lmodels_tokeep]
  lweight_tokeep <- CvM_df_filtered$weight
  
  p12_multimodel <- multimodel_average(lp12_tokeep, lmodels_tokeep, lweight = lweight_tokeep)  
  
  list(
    lp12 = lp12,
    lUhat = lUhat,
    CvM_df = CvM_df,
    p12_multimodel = p12_multimodel
  )
}
  
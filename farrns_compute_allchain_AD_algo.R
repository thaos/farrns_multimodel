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
  lmodels <- lmodels[lmodels != "HadCRUT"]
  lUhat <- lapply(lmodels, function(model, tas_cmip5, tas_hadcrut, kernel, bandwidth){
    # browser()
    print("---------------------------")
    print(model)
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" & tas_model$year <= 1900, ]
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
    Uhat <- compute_Uhat(
      Xm = tas_model_counterfactual$tas,
      Zm = tas_model_factual$tas
    )
    return(Uhat)
  }, tas_cmip5 = tas_cmip5, tas_hadcrut = tas_hadcrut, kernel = kernel, bandwidth = bandwidth)
  
  CvM_df <- CvM_checkH0(
    lUhat = lUhat,
    lmodels = lmodels
  )
  CvM_df_filtered <- keep_onemodel_perinstitute(CvM_df, unique(tas_cmip5[, c("institute", "model")]))
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
  

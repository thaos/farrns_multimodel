compute_allchain <- function(lmodels, linstitutes, tas_cmip5, kernel, bandwidth) {
  lp12 <- lapply(lmodels, function(model, tas_cmip5, kernel, bandwidth) {
    # print("---------------------------")
    # print(model)
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" | tas_model$experiment == "rcp85", ]
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
    # print(head(tas_model_counterfactual))
    # tas_model_counterfactual <- tas_model[tas_model$experiment == "historical" & tas_model$year <= 1900, ]
    tpred <- min(tas_model$year):max(tas_model$year)
    p12 <- estim_p12_ns(x = tas_model_counterfactual$tas, t = tas_model_factual$year, z = tas_model_factual$tas, tpred = tpred, kernel = kernel, h = bandwidth)
  }, tas_cmip5 = tas_cmip5, kernel = kernel, bandwidth = bandwidth)

  lUhat <- lapply(lmodels, function(model, tas_cmip5, tas_hadcrut, kernel, bandwidth) {
    # browser()
    # print("---------------------------")
    # print(model)
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
  CvM_df <- keep_onemodel_perinstitute(
    CvM_df,
    linstitutes =  linstitutes
  )

  imodels_tokeep <- which(CvM_df$weight > 0)
  lmodels_tokeep <- CvM_df$model[imodels_tokeep]
  lp12_tokeep <- lp12[imodels_tokeep]
  expert_weights <- compute_expert_weights(lp12_tokeep, lmodels_tokeep)
  lweight_tokeep <- expert_weights$weight
  CvM_df$weight[imodels_tokeep] <- lweight_tokeep

  p12_multimodel <- multimodel_average(lp12_tokeep, lmodels_tokeep, lweight = lweight_tokeep)

  list(
    lp12 = lp12,
    lUhat = lUhat,
    CvM_df = CvM_df,
    p12_multimodel = p12_multimodel
  )
}
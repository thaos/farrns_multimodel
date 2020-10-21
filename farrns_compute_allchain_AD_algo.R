compute_allchain <- function(lmodels, linstitutes, tas_cmip5, kernel, bandwidth) {
  lp12 <- lapply(lmodels, function(model, tas_cmip5, kernel, bandwidth) {
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" | tas_model$experiment == "rcp85", ]
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
    tpred <- min(tas_model$year):max(tas_model$year)
    p12 <- estim_p12_ns(x = tas_model_counterfactual$tas, t = tas_model_factual$year, z = tas_model_factual$tas, tpred = tpred, kernel = kernel, h = bandwidth)
  }, tas_cmip5 = tas_cmip5, kernel = kernel, bandwidth = bandwidth)

  lXmZm <- lapply(lmodels, function(model, tas_cmip5) {
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" & tas_model$year <= 1900, ]
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat" & tas_model$year <= 1900, ]
    return(
      list(
        Xm = tas_model_counterfactual$tas,
        Zm = tas_model_factual$tas
      )
    )
  }, tas_cmip5 = tas_cmip5)

  lweights <- compute_expert_weights(lp12)

  p12_multimodel <- multimodel_average(lp12, lmodels, lweights = lweights)

  list(
    lp12 = lp12,
    lXmZm = lXmZm,
    lweights = lweights,
    p12_multimodel = p12_multimodel
  )
}

reformat_lp12 <- function(lp12_pergridpoint, lmodels, linstitutes) {
  p12 <- mapply(function(p12_onegridpoint, igridpoint) {
    p12_multimodel_df <- data.frame(
      igridpoint = igridpoint,
      institute = "multimodel",
      model = "multimodel",
      year = p12_onegridpoint$p12_multimodel$year,
      p12 = p12_onegridpoint$p12_multimodel$p12,
      sig = p12_onegridpoint$p12_multimodel$sig
    )
    lp12_model_df <- mapply(function(p12_df, model, institute) {
      p12_model_df <- data.frame(
        igridpoint = igridpoint,
        institute = institute,
        model = model,
        year = p12_df$tpred,
        p12 = p12_df$p12_hat,
        sig = p12_df$sigma_p12_hat
      )
    },
    p12_df = p12_onegridpoint$lp12, model = lmodels, institute = linstitutes,
    SIMPLIFY = FALSE
    ) %>% do.call(rbind, .)
    rbind(lp12_model_df, p12_multimodel_df)
  },
  p12_onegridpoint = lp12_pergridpoint,
  igridpoint = iworld,
  SIMPLIFY = FALSE
  ) %>% 
  do.call(rbind, .)
}

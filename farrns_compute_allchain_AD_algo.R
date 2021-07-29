compute_allchain <- function(lmodels, linstitutes, tas_cmip5, kernel, bandwidth,year_pi = c(1850, 1900)) {
  lp12 <- lapply(lmodels, function(model, tas_cmip5, kernel, bandwidth) {
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" | tas_model$experiment == "rcp85", ]
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
    tpred <- min(tas_model$year):max(tas_model$year)
    p12 <- estim_p12_ns(x = tas_model_counterfactual$tas, t = tas_model_factual$year, z = tas_model_factual$tas, tpred = tpred, kernel = kernel, h = bandwidth)
  }, tas_cmip5 = tas_cmip5, kernel = kernel, bandwidth = bandwidth)

  lXmZm <- lapply(lmodels, function(model, tas_cmip5) {
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" & tas_model$year >= year_pi[1] & tas_model$year <= year_pi[2], ]
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat" & tas_model$year >= year_pi[1] & tas_model$year <= year_pi[2], ]
    return(
      list(
        Xm = tas_model_counterfactual$tas,
        Zm = tas_model_factual$tas
      )
    )
  }, tas_cmip5 = tas_cmip5)
  
  pvalue <- sapply(lXmZm, function(XmZm) twosamples::ad_test(XmZm$Xm, XmZm$Zm)[2])
  imodels <- pvalue > 0.2
  
  dates_H0 <- seq.int(year_pi[1], year_pi[2])
  experts <- vapply(
    lp12,
    function(p12) {
      return(
        p12$p12_hat[p12$tpred >= min(dates_H0) & p12$tpred <= max(dates_H0)]
      )
    },
    FUN.VALUE = numeric(length(dates_H0))
  )

  kl <- apply(experts, 2, kldiv)
  lweights_kl <- rep(0, length(lp12))
  lweights_oracle <- rep(0, length(lp12))
  lweights_kl[imodels] <- compute_kl_weights(experts[, imodels])
  lweights_oracle[imodels] <- compute_oracle_weights(experts[, imodels])
  p12_multimodel_kl <- multimodel_average(lp12, lmodels, lweights = lweights_kl)
  p12_multimodel_oracle <- multimodel_average(lp12, lmodels, lweights = lweights_oracle)
  lweights_best <- kl == min(kl)
  p12_multimodel_best <- multimodel_average(lp12, lmodels, lweights = lweights_best)

  list(
    lp12 = lp12,
    lXmZm = lXmZm,
    lweights_kl = lweights_kl,
    lweights_oracle = lweights_oracle,
    lweights_best = lweights_best,
    p12_multimodel_kl = p12_multimodel_kl,
    p12_multimodel_oracle = p12_multimodel_oracle,
    p12_multimodel_best = p12_multimodel_best
  )
}

reformat_lp12 <- function(lp12_pergridpoint, lmodels, linstitutes) {
  p12 <- mapply(function(p12_onegridpoint, igridpoint) {
    p12_multimodel_kl_df <- data.frame(
      igridpoint = igridpoint,
      institute = "multimodel_kl",
      model = "multimodel_kl",
      year = p12_onegridpoint$p12_multimodel_kl$year,
      p12 = p12_onegridpoint$p12_multimodel_kl$p12,
      sig = p12_onegridpoint$p12_multimodel_kl$sig
    )
    p12_multimodel_oracle_df <- data.frame(
      igridpoint = igridpoint,
      institute = "multimodel_oracle",
      model = "multimodel_oracle",
      year = p12_onegridpoint$p12_multimodel_oracle$year,
      p12 = p12_onegridpoint$p12_multimodel_oracle$p12,
      sig = p12_onegridpoint$p12_multimodel_oracle$sig
    )
    p12_multimodel_best_df <- data.frame(
      igridpoint = igridpoint,
      institute = "multimodel_best",
      model = "multimodel_best",
      year = p12_onegridpoint$p12_multimodel_best$year,
      p12 = p12_onegridpoint$p12_multimodel_best$p12,
      sig = p12_onegridpoint$p12_multimodel_best$sig
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
    rbind(lp12_model_df, p12_multimodel_kl_df, p12_multimodel_oracle_df, p12_multimodel_best_df)
  },
  p12_onegridpoint = lp12_pergridpoint,
  igridpoint = iworld,
  SIMPLIFY = FALSE
  ) %>% 
  do.call(rbind, .)
}

library(magrittr)
library(twosamples)

ecdf_nw <- function(q, tpred, t, x, kernel, bandwidth){
  dvect <- abs(tpred - t)
  kvect <- kernel(dvect, h = bandwidth)
  indicator <- (x < q)
  weighted.mean(indicator, w = kvect)
}

quantile_nw <- function(p, tpred, t, x, kernel, bandwidth){
  # cat("p = ", p , "; tpred = ", tpred,  " \n")
  dvect <- abs(tpred - t)
  kvect <- kernel(dvect, h = bandwidth)
  uniroot(f = function(q) {
    indicator <- (x < q)
    Gq <- weighted.mean(indicator, w = kvect)
    Gq - p
  }, interval = range(x) + c(-sd(x), +sd(x)))$root
}

compute_Uhat <- function(t, tm, Xm, Z, Zm, kernel, bandwidth){
  FmXm <- sapply(seq_along(tm), function(i){
    ecdf_nw(Xm[i], tm[i], tm, Zm, kernel = kernel_epanechnikov, bandwidth = bandwidth)
  })
  Uhat <- sapply(seq_along(tm), function(i){
    quantile_nw(FmXm[i], tm[i], t, Z, kernel = kernel_epanechnikov, bandwidth = bandwidth)
  })
  return(Uhat)
}

check_HO_qqplot <- function(X, Uhat){
  xyrange <- range(X, Uhat)
  qqplot(X, Uhat, main = "", xlim = xyrange, ylim = xyrange)
  abline(a = 0, b = 1, col = "red")
}

check_HO_rankhist <- function(X, Uhat){
  GmU <- ecdf(X)(Uhat)
  cut_GmU <- cut(GmU, seq(0, 1, 0.2), include.lowest = TRUE)
  # ks.test(GmU[GmU < 1 & GmU > 0], y = "punif")
  print(chisq.test(table(cut_GmU)))
  hist(GmU, breaks = 5, main = "")
  abline(h = length(GmU)/5, col = "red")
}

hist_checkH0 <- function(X, lUhat, lmodels, linstitutes){
  GmU_df <- mapply(
    function(Uhat, model, institute){
      GmU <- ecdf(X)(Uhat)
      data.frame(institute = institute, model = model, GmU = GmU)
    },
    Uhat = lUhat, model = lmodels, institute = linstitutes,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  print(head(GmU_df))
  ggplot(data = GmU_df, aes(x = GmU, fill = institute, col = grepl("(CM6)", model))) + 
    geom_histogram(aes(y = ..density..), breaks = seq(0, 1, by = 0.2)) + 
    geom_hline(yintercept = 1) +
    facet_wrap( ~ institute + model, ncol = 5) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("white", "black"), labels = c("CMIP5", "CMIP6"), name= "CMIP")+
    ggtitle("Checking H0: ranks of Uhat with respect to X")
}

qqplot_checkH0 <- function(X, lUhat, lmodels, linstitutes){
  qq_df <- mapply(
    function(Uhat, model, institute){
      qqdata <- qqplot(X, Uhat, plot.it = FALSE) %>% as.data.frame()
      names(qqdata) <- c("X", "Uhat")
      cbind(institute = institute, model = model, qqdata)
    },
    Uhat = lUhat, model = lmodels, institute = linstitutes,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  print(head(qq_df))
  xylim = range(qq_df$X, qq_df$Uhat)
  ggplot(data = qq_df) + 
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x = X, y = Uhat, col = institute)) +
    facet_wrap( ~ institute + model, ncol = 5) +
    theme(legend.position = "none") +
    ggtitle("Checking H0: qq-plot(X, Uhat)") +
    coord_fixed(xlim = xylim, ylim = xylim) 
}

CvM_checkH0 <- function(X, lUhat, lmodels){
  CvM_df <- mapply(
    function(Uhat, model){
      CvM <- CramerVonMisesTwoSamples(X, Uhat)
      data.frame(model = model, CvM = CvM)
    },
    Uhat = lUhat, model = lmodels,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  # CvM_df$CvM[lmodels != "HadCRUT"] <- pmax(CvM_df$CvM[lmodels != "HadCRUT"] - CvM_df$CvM[lmodels == "HadCRUT"],  CvM_df$CvM[lmodels == "HadCRUT"])
  # CvM_df$CvM[lmodels != "HadCRUT"] <- CvM_df$CvM[lmodels != "HadCRUT"] * CvM_df$CvM[lmodels == "HadCRUT"] /  min(CvM_df$CvM[lmodels != "HadCRUT"])
  CvM_df$CvM[lmodels == "HadCRUT"] <- min(CvM_df$CvM[lmodels != "HadCRUT"])
  CvM_df$weight <- (1/CvM_df$CvM) / sum(1/CvM_df$CvM) 
  return(CvM_df)
}


keep_onemodel_perinstitute <- function(CvM_df, InstituteModel_df){
  CvM_df <- merge(x = InstituteModel_df, y = CvM_df, by = "model")
  CvM_df <- split.data.frame(CvM_df, f = CvM_df$institute) %>% 
    lapply(., function(df){
      print(df$institute)
      if(nrow(df)> 1){
       df[(min(df$CvM) != df$CvM) , "weight"] <- 0
      }
      return(df)
    }) %>%
    do.call(rbind, .)
  return(CvM_df)
}

multimodel_average <- function(lp12, lmodels, lweight = 1/length(lmodels)){
  p12_df <- mapply(function(p12, model, weight){
    data.frame( 
      model = model,
      weight = weight, 
      year = p12$tpred,
      p12 = p12$p12_hat,
      sig = p12$sigma_p12_hat
    )
  }, p12 = lp12, model = lmodels, weight = lweight, SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  print(head(p12_df))
  p12_byyear <- split(p12_df, f = p12_df$year)
  lapply(p12_byyear, function(df){
    weight <- df$weight / sum(df$weight)
    p12 <- sum(df$p12 * weight)
    var_intra <- sum(df$sig^2 * weight^2)
    var_inter <- sum((df$p12 - p12)^2 * weight^2)
    print(var_inter/(var_intra + var_inter))
    sig <- sqrt(var_intra + var_inter)
    data.frame(year = unique(df$year), p12  = p12, sig = sig)
  }) %>% do.call(rbind, .)
}

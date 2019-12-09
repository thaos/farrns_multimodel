library(magrittr)
library(twosamples)
library(goftest)


compute_Uhat <- function(Xm, Zm){
  Uhat <- ecdf(Xm)(Zm)
}

check_HO_qqplot <- function(X, Uhat){
  xyrange <- range(X, Uhat)
  qqplot(X, Uhat, main = "", xlim = xyrange, ylim = xyrange)
  abline(a = 0, b = 1, col = "red")
}

check_HO_rankhist <- function(Uhat){
  cut_Uhat <- cut(Uhat, seq(0, 1, 0.2), include.lowest = TRUE)
  # ks.test(GmU[GmU < 1 & GmU > 0], y = "punif")
  print(chisq.test(table(cut_Uhat)))
  hist(Uhat, breaks = 5, main = "")
  abline(h = length(Uhat)/5, col = "red")
}

hist_checkH0 <- function(X, lUhat, lmodels, linstitutes){
  Uhat_df <- mapply(
    function(Uhat, model, institute){
      data.frame(institute = institute, model = model, Uhat = Uhat)
    },
    Uhat = lUhat, model = lmodels, institute = linstitutes,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  ggplot(data = Uhat_df, aes(x = Uhat, fill = institute, col = grepl("(CM6)", model))) + 
    geom_histogram(aes(y = ..density..), breaks = seq(0, 1, by = 0.2)) + 
    geom_hline(yintercept = 1) +
    facet_wrap( ~ institute + model, ncol = 5) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("white", "black"), labels = c("CMIP5", "CMIP6"), name= "CMIP")+
    ggtitle("Checking H0: Uhat histogram")
}

qqplot_checkH0 <- function(X, lUhat, lmodels, linstitutes){
  qq_df <- mapply(
    function(Uhat, model, institute){
      qqdata <- data.frame(
        qunif = qunif(p = rank(Uhat)/length(Uhat)),
        Uhat = Uhat
      )
      cbind(institute = institute, model = model, qqdata)
    },
    Uhat = lUhat, model = lmodels, institute = linstitutes,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  print(head(qq_df))
  xylim = range(qq_df$qunif, qq_df$Uhat)
  ggplot(data = qq_df) + 
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x = qunif, y = Uhat, col = institute)) +
    facet_wrap( ~ institute + model, ncol = 5) +
    theme(legend.position = "none") +
    ggtitle("Checking H0: qq-plot(qunif, Uhat)") +
    coord_fixed(xlim = xylim, ylim = xylim) 
}

CvM_checkH0 <- function(lUhat, lmodels){
  CvM_df <- mapply(
    function(Uhat, model){
      CvM <- cvm.test(Uhat, "punif")
      data.frame(model = model, CvM = CvM$statistic)
    },
    Uhat = lUhat, model = lmodels,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
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

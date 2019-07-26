library(ncdf4)
library(magrittr)
library(ggplot2)
library(SpecsVerification)
library(twosamples)
devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("~/2_Code/Naveau/farFAR/farr/inprogress/tests_nsfar_algo.R")
source("check_H0_algo.R")

ihadcrut <- length(lmodels)

bandwidth <- theta$h
par(mfrow = c(1, 1), cex.main = 0.6, cex.lab = 0.6)
pdf("check_H0_allruns.pdf")
lUhat <- lapply(lmodels, function(model, tas_cmip5, tas_hadcrut, kernel, bandwidth){
  # browser()
  print("---------------------------")
  print(model)
  tas_model <- tas_cmip5[tas_cmip5$model == model, ]
  tas_model_factual <- tas_model[tas_model$experiment == "historical", ]
  if(model == "HadCRUT"){
    tas_model_counterfactual <- tas_model[tas_model$year <= 1900 , ]
  } else {
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
  }
  yend <- min(c(max(tas_hadcrut$year), max(tas_model_factual$year), max(tas_model_counterfactual$year)))
  ystart <- max(c(min(tas_hadcrut$year), min(tas_model_factual$year), min(tas_model_counterfactual$year)))
  tas_hadcrut <- tas_hadcrut[tas_hadcrut$year <= yend & tas_hadcrut$year >= ystart, ]
  tas_model_counterfactual <- tas_model_counterfactual[tas_model_counterfactual$year <= yend  & tas_model_counterfactual$year >= ystart, ]
  tas_model_factual <- tas_model_factual[tas_model_factual$year <= yend  & tas_model_factual$year >= ystart, ]
  Uhat <- compute_Uhat(
    t = tas_hadcrut$year,
    tm = tas_model_counterfactual$year,
    Xm = tas_model_counterfactual$tas,
    Z = tas_hadcrut$tas,
    Zm = tas_model_factual$tas,
    kernel = kernel, bandwidth = bandwidth
  )
  check_HO_rankhist(tas_hadcrut_counterfactual$tas, Uhat)
  title(main = model)
  check_HO_qqplot(tas_hadcrut_counterfactual$tas, Uhat)
  title(main = model)
  return(Uhat)
}, tas_cmip5 = tas_cmip5, tas_hadcrut = tas_hadcrut, kernel = krnl, bandwidth = bandwidth)
dev.off()

hist_checkH0 <- function(X, lUhat, lmodels){
  GmU_df <- mapply(
    function(Uhat, model){
      GmU <- ecdf(X)(Uhat)
      data.frame(model = model, GmU = GmU)
    },
    Uhat = lUhat, model = lmodels,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  print(head(GmU_df))
  ggplot(data = GmU_df, aes(x = GmU, fill = model)) + 
    geom_histogram(aes(y = ..density..), breaks = seq(0, 1, by = 0.2)) + 
    geom_hline(yintercept = 1) +
    facet_wrap( ~ model, ncol = 5) +
    theme(legend.position = "none") +
    ggtitle("Checking H0: ranks of Uhat with respect to X")
}
hist_checkH0( 
  X = tas_hadcrut_counterfactual$tas,
  lUhat = lUhat,
  lmodels = lmodels[-ihadcrut]
)
ggsave("hist_checkH0_ggplot.pdf", dpi = "retina", width = 20, height = 15, units = "cm")

qqplot_checkH0 <- function(X, lUhat, lmodels){
  qq_df <- mapply(
    function(Uhat, model){
      qqdata <- qqplot(X, Uhat, plot.it = FALSE) %>% as.data.frame()
      names(qqdata) <- c("X", "Uhat")
      cbind(model = model, qqdata)
    },
    Uhat = lUhat, model = lmodels,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  print(head(qq_df))
  ggplot(data = qq_df) + 
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x = X, y = Uhat, col = model)) +
    facet_wrap( ~ model, ncol = 5) +
    theme(legend.position = "none") +
    ggtitle("Checking H0: qq-plot(X, Uhat)")
}
qqplot_checkH0(
  X = tas_hadcrut_counterfactual$tas,
  lUhat = lUhat,
  lmodels = lmodels
)
ggsave("qqplot_checkH0_ggplot.pdf", dpi = "retina", width = 20, height = 15, units = "cm")

CvM_checkH0 <- function(X, lUhat, lmodels){
  CvM_df <- mapply(
    function(Uhat, model){
      CvM <- cvm_stat(X, Uhat)
      data.frame(model = model, CvM = CvM)
    },
    Uhat = lUhat, model = lmodels,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  CvM_df$CvM[lmodels != "HadCRUT"] <- pmax(CvM_df$CvM[lmodels != "HadCRUT"] - CvM_df$CvM[lmodels == "HadCRUT"],  CvM_df$CvM[lmodels == "HadCRUT"])
  CvM_df$CvM[lmodels != "HadCRUT"] <- CvM_df$CvM[lmodels != "HadCRUT"] * CvM_df$CvM[lmodels == "HadCRUT"] /  min(CvM_df$CvM[lmodels != "HadCRUT"])
  CvM_df$weight <- (1/CvM_df$CvM) / sum(1/CvM_df$CvM) 
  return(CvM_df)
}
CvM_df <- CvM_checkH0(
  X = tas_hadcrut_counterfactual$tas,
  lUhat = lUhat,
  lmodels = lmodels
)

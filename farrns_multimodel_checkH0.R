library(ncdf4)
library(magrittr)
library(ggplot2)
devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("~/2_Code/Naveau/farFAR/farr/inprogress/tests_nsfar_algo.R")
source("check_H0_algo.R")

tas_cmip5 <- readRDS(file = "tas_cmip5.rds")
nbtimestep <- aggregate(tas ~ model + experiment, data = tas_cmip5, FUN = length)

lmodels <- unique(tas_cmip5$model) %>% as.character()
lmodels <- subset(nbtimestep, experiment == "rcp85", select = model) %>% unlist() %>% as.character()

tas_cmip5 <- tas_cmip5[tas_cmip5$model %in% lmodels, ]
# tas_cmip5 <- subset(tas_cmip5, run == "r1i1p1")

nc = nc_open(file = "tas_hadcrut.nc")
tas = ncvar_get(nc, "temperature_anomaly") %>% as.numeric()
year = ncvar_get(nc, "time") %/% 10000 %>% as.numeric()
tas_hadcrut <- data.frame("institute" = "HadCRUT", "model" = "HadCRUT", "experiment" = "historical", "run" = "obs", "year" = year,  "tas" = tas)
# 
# tas_cmip5 <- rbind(tas_hadcrut, tas_cmip5)
lmodels <- unique(tas_cmip5$model) %>% as.character()

krnl <- kernel_epanechnikov
# knrl <- kernel_gauss


tas_hadcrut_counterfactual <- subset(tas_hadcrut, year <= 1900)

# first fit to get bandwidth h
theta <- estim_theta.nswexp(x = tas_hadcrut_counterfactual$tas, t = tas_hadcrut$year, z = tas_hadcrut$tas, kernel = krnl, h = NULL)
print(theta$utest_pvalue)
par(mfrow = c(1, 3), cex.main = 0.6, cex.lab = 0.6)
hist(theta) 
ecdf(theta)
qqplot(theta)

bandwidth <- theta$h
par(mfrow = c(1, 1), cex.main = 0.6, cex.lab = 0.6)
pdf("check_H0_allruns.pdf")
lUhat <- lapply(lmodels, function(model, tas_cmip5, tas_hadcrut, kernel, bandwidth){
  # browser()
  print("---------------------------")
  print(model)
  tas_hadcrut_counterfactual <- subset(tas_hadcrut, year <= 1900)
  tas_model <- tas_cmip5[tas_cmip5$model == model, ]
  tas_model_factual <- tas_model[tas_model$experiment == "historical", ]
  tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
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

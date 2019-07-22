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
tas_cmip5 <- rbind(tas_hadcrut, tas_cmip5)
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
lmodels_tokeep <- c(
  "bcc-csm1-1", "CanESM2","CNRM-CM5", "CSIRO-Mk3-6-0",
  "IPSL-CM5A-LR", "MRI-CGCM3",  "GISS-E2-H", "GISS-E2-R",
  "CCSM4",  "NorESM1-M", "CESM1-CAM5", "HadCRUT"
)

gen_boot_pool <- function(lmodels, tas_cmip5, bootstrap = FALSE){
  tas_cmip5_factual <- tas_cmip5[tas_cmip5$experiment != "historicalNat" & tas_cmip5$year <= 2100, ]
  tas_cmip5_counterfactual <- tas_cmip5[tas_cmip5$experiment == "historicalNat", ]
  if(bootstrap == TRUE){
    tas_cmip5_counterfactual <- tas_cmip5_counterfactual[sample.int(nrow(tas_cmip5_counterfactual), replace = TRUE), ]
    tas_cmip5_factual <- tas_cmip5_factual[sample.int(nrow(tas_cmip5_factual), replace = TRUE), ]
  }
  lGmZ <- lapply(lmodels, function(model, tas_cmip5_factual, tas_cmip5_counterfactual){
    # browser()
    print("---------------------------")
    print(model)
    if(model == "HadCrut"){
      x <- tas_cmip5_factual$tas[tas_cmip5_factual$year <= 1900]
    } else {
      x = tas_cmip5_counterfactual$tas
    }
    z = tas_cmip5_factual$tas
    t = tas_cmip5_factual$year
    Gm <- ecdf(x)
    GmZ <- Gm(z)
    data.frame(year = t, GmZ = GmZ)
  }, tas_cmip5_factual = tas_cmip5_factual, tas_cmip5_counterfactual = tas_cmip5_counterfactual)
  GmZ_df <- do.call(rbind, lGmZ)
}
boot_GmZ_df <- lapply(seq.int(250), function(i, lmodels, tas_cmip5) {
  if (i == 1) gen_boot_pool(lmodels, tas_cmip5, bootstrap = FALSE)
  else gen_boot_pool(lmodels, tas_cmip5, bootstrap = TRUE)
}, lmodels = lmodels_tokeep, tas_cmip5 = tas_cmip5)

estim_p12_from_GmZ <- function(tpred, t, GmZ, kernel, bandwidth){
  p12 <- sapply(tpred, function(ti){
    dvect <- abs(ti - t)
    kvect <- kernel(dvect, h = bandwidth)
    weighted.mean(GmZ, w = kvect)
  })
  return(p12)
}
p12_pool_boot <- sapply(boot_GmZ_df, function(df, tpred){
  estim_p12_from_GmZ(tpred = tpred, t = df$year, GmZ = df$GmZ, kernel = kernel_epanechnikov, bandwidth = bandwidth) 
}, tpred = 1850:2100)
ci95_pool <- apply(p12_pool_boot, 1, quantile, probs = c(0.05, 0.5, 0.95))
par(mfrow = c(1,1))
plot(boot_GmZ_df[[1]]$year, boot_GmZ_df[[1]]$GmZ, col = "grey",  cex = 0.1)
matlines(1850:2100, t(ci95_pool), lty = c(2, 1, 2), lwd = c(1, 2, 1), col = "black")

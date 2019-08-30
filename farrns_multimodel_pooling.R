library(ncdf4)
library(magrittr)
library(ggplot2)
devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("~/2_Code/Naveau/farFAR/farr/inprogress/tests_nsfar_algo.R")
source("check_H0_algo.R")

bandwidth <- theta$h
lmodels_tokeep <- c(
  "bcc-csm1-1", "CanESM2","CNRM-CM5", "CSIRO-Mk3-6-0",
  "IPSL-CM5A-LR", "MRI-CGCM3",  "GISS-E2-R", 
  "CCSM4",  "NorESM1-M", "CESM1-CAM5"#, "HadCRUT"
)
lmodels_tokeep <- c(
  "bcc-csm1-1", "CanESM2","CNRM-CM5", "CSIRO-Mk3-6-0",
  "IPSL-CM5A-LR", "MIROC-ESM-CHEM", "MRI-CGCM3",  "GISS-E2-R",
  "CCSM4",  "NorESM1-M", "GFDL-CM3", "CESM1-CAM5", "HadCRUT"
)
# lmodels_tokeep <- lmodels

lp12_tokeep <- lp12[lmodels %in% lmodels_tokeep]
lweight_tokeep <- CvM_df$weight[lmodels %in% lmodels_tokeep]

CvM_df_filtered <- keep_onemodel_perinstitute(CvM_df, unique(tas_cmip5[, c("institute", "model")]))
lmodels_tokeep <- CvM_df_filtered$model
lp12_tokeep <- lp12[lmodels %in% lmodels_tokeep]
lweight_tokeep <- CvM_df_filtered$weight

p12_multimodel_ggplot <- multimodel_average(lp12_tokeep, lmodels_tokeep)  
p12_multimodel_ggplot <- multimodel_average(lp12_tokeep, lmodels_tokeep, lweight = lweight_tokeep)  
p <- ggplot(p12_multimodel_ggplot) +
  geom_hline(yintercept = 1/2) +
  geom_line(aes(y = p12, x = year), lwd = 1.2) +
  geom_ribbon(aes(ymin = p12 - qnorm(0.95) * sig, ymax = p12 + qnorm(0.95) * sig, x = year), alpha = 0.3) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) + ggtitle("q1(t), multimodel synthesis") + ylab("q1")
plot(p)
ggsave("p12_multimodel_cvm_ggplot.pdf", dpi = "retina", width = 20, height = 15, units = "cm")

p12_allinone_ggplot <- rbind(
  lp12boot_ggplot,
  cbind(
    model = "Multi-Model",
    p12_multimodel_ggplot[, 1:2],
    ci_q05 = p12_multimodel_ggplot[, 2] - qnorm(0.95) * p12_multimodel_ggplot[, 3],
    ci_q95 = p12_multimodel_ggplot[, 2] + qnorm(0.95) * p12_multimodel_ggplot[, 3]
  )
)
p <- ggplot(p12_allinone_ggplot) +
  geom_hline(yintercept = 1/2) +
  geom_line(aes(y = p12, x = year, colour = model), lwd = 1.2) +
  geom_ribbon(aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = model), alpha = 0.3) +
  facet_wrap( ~ model, ncol = 5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +  ylim(min(lp12boot_ggplot$ci_q95), 1) + 
  ggtitle("q1(t), factual = historical + rcp85, counterfactual = historicalNat") + ylab("q1")
plot(p)
ggsave("p12_allinone_cvm_ggplot.pdf", dpi = "retina", width = 20, height = 15, units = "cm")

# gen_boot_pool <- function(lmodels, tas_cmip5, bootstrap = FALSE){
#   tas_cmip5_factual <- tas_cmip5[tas_cmip5$experiment != "historicalNat" & tas_cmip5$year <= 2100, ]
#   tas_cmip5_counterfactual <- tas_cmip5[tas_cmip5$experiment == "historicalNat", ]
#   if(bootstrap == TRUE){
#     tas_cmip5_counterfactual <- tas_cmip5_counterfactual[sample.int(nrow(tas_cmip5_counterfactual), replace = TRUE), ]
#     tas_cmip5_factual <- tas_cmip5_factual[sample.int(nrow(tas_cmip5_factual), replace = TRUE), ]
#   }
#   lGmZ <- lapply(lmodels, function(model, tas_cmip5_factual, tas_cmip5_counterfactual){
#     # browser()
#     print("---------------------------")
#     print(model)
#     if(model == "HadCrut"){
#       x <- tas_cmip5_factual$tas[tas_cmip5_factual$year <= 1900]
#     } else {
#       x = tas_cmip5_counterfactual$tas
#     }
#     z = tas_cmip5_factual$tas
#     t = tas_cmip5_factual$year
#     Gm <- ecdf(x)
#     GmZ <- Gm(z)
#     data.frame(year = t, GmZ = GmZ)
#   }, tas_cmip5_factual = tas_cmip5_factual, tas_cmip5_counterfactual = tas_cmip5_counterfactual)
#   GmZ_df <- do.call(rbind, lGmZ)
# }
# boot_GmZ_df <- lapply(seq.int(250), function(i, lmodels, tas_cmip5) {
#   if (i == 1) gen_boot_pool(lmodels, tas_cmip5, bootstrap = FALSE)
#   else gen_boot_pool(lmodels, tas_cmip5, bootstrap = TRUE)
# }, lmodels = lmodels_tokeep, tas_cmip5 = tas_cmip5)
# 
# estim_p12_from_GmZ <- function(tpred, t, GmZ, kernel, bandwidth){
#   p12 <- sapply(tpred, function(ti){
#     dvect <- abs(ti - t)
#     kvect <- kernel(dvect, h = bandwidth)
#     weighted.mean(GmZ, w = kvect)
#   })
#   return(p12)
# }
# p12_pool_boot <- sapply(boot_GmZ_df, function(df, tpred){
#   estim_p12_from_GmZ(tpred = tpred, t = df$year, GmZ = df$GmZ, kernel = kernel_epanechnikov, bandwidth = bandwidth) 
# }, tpred = 1850:2100)
# ci95_pool <- apply(p12_pool_boot, 1, quantile, probs = c(0.05, 0.5, 0.95))
# par(mfrow = c(1,1))
# plot(boot_GmZ_df[[1]]$year, boot_GmZ_df[[1]]$GmZ, col = "grey",  cex = 0.1)
# matlines(1850:2100, t(ci95_pool), lty = c(2, 1, 2), lwd = c(1, 2, 1), col = "black")
# abline(h = 0.5, col = "red")

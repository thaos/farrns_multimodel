library(ncdf4)
library(magrittr)
library(ggplot2)
devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("~/2_Code/Naveau/farFAR/farr/inprogress/tests_nsfar_algo.R")

tas_cmip5 <- readRDS(file = "tas_cmip5.rds")
nbtimestep <- aggregate(tas ~ model + experiment, data = tas_cmip5, FUN = length)

lmodels <- unique(tas_cmip5$model) %>% as.character()
lmodels <- subset(nbtimestep, experiment == "rcp85", select = model) %>% unlist() %>% as.character()

tas_cmip5 <- tas_cmip5[tas_cmip5$model %in% lmodels, ]
tas_cnrm <- tas_cmip5[tas_cmip5$model == "CNRM-CM5", ]
tas_cnrm_factual <- tas_cmip5[tas_cmip5$experiment != "historicalNat", ]
tas_cnrm_counterfactual <- tas_cmip5[tas_cmip5$experiment == "historicalNat", ]
# tas_cmip5 <- subset(tas_cmip5, run == "r1i1p1")

write.table(
  aggregate(
    run ~ institute + model, 
    FUN = function(x) paste(sort(unique(x)), collapse = ", "),
    data = tas_cmip5
  ),
  file = "list_models_runs.txt",
  row.names = FALSE
)


nc = nc_open(file = "tas_hadcrut.nc")
tas = ncvar_get(nc, "temperature_anomaly") %>% as.numeric()
year = ncvar_get(nc, "time") %/% 10000 %>% as.numeric()
nc_close(nc)
tas_hadcrut <- data.frame("institute" = "HadCRUT", "model" = "HadCRUT", "experiment" = "historical", "run" = "obs", "year" = year,  "tas" = tas)
# 
tas_cmip5 <- rbind(tas_cmip5, tas_hadcrut)
lmodels <- unique(tas_cmip5$model) %>% as.character()

krnl <- kernel_epanechnikov
# knrl <- kernel_gauss


tas_hadcrut_counterfactual <- subset(tas_hadcrut, year <= 1900)

# first fit to get bandwidth h
theta <- estim_theta.nswexp(x = tas_hadcrut_counterfactual$tas, t = tas_hadcrut$year, z = tas_hadcrut$tas, kernel = krnl, h = NULL)
# theta <- estim_theta.nswexp(x = tas_cnrm_counterfactual$tas, t = tas_cnrm_factual$year, z = tas_cnrm_factual$tas, kernel = krnl, h = NULL)

bandwidth <- theta$h
# bandwidth <- 90

thetaboot <- boot_theta_fit.nswexp(x = tas_hadcrut_counterfactual$tas, t = tas_hadcrut$year, z = tas_hadcrut$tas, kernel = krnl, h = bandwidth, B = 250)
p12boot <- boot_p1r_fit.nswexp(thetaboot, rp = 2)
p12_boot_mean <- apply(p12boot$p1r_boot[, 1,], 1, mean)
p12_boot_ci90 <- t(apply(p12boot$p1r_boot[, 1,], 1, quantile, probs = c(0.05, 0.95)))
p12_ci90 <- data.frame(year = thetaboot$t_unique, q05 = p12_boot_ci90[, 1], mean = p12_boot_mean, q95 = p12_boot_ci90[, 2]) %>% unique()
p12_ci90[order(p12_ci90$year), ]

# Comparison of Asymptotic and Bootstrap CI
p12 <- estim_p12_ns(x = tas_hadcrut_counterfactual$tas, t = tas_hadcrut$year, z = tas_hadcrut$tas, kernel = krnl, h = bandwidth)
with(p12, plot(t_unique, p12_hat, lty = 2, ylim = range(p12_ci90[, -1], p12_hat + 1.64 * sigma_p12_hat, p12_hat - 1.64 * sigma_p12_hat)))
with(p12, lines(t_unique, p12_hat + 1.64 * sigma_p12_hat))
with(p12, lines(t_unique, p12_hat - 1.64 * sigma_p12_hat))
with(p12_ci90, points(year, mean, col = "green"))
with(p12_ci90, lines(year, q05, col = "green"))
with(p12_ci90, lines(year, q95, col = "green"))

lp12 <- lapply(lmodels, function(model, tas_cmip5, kernel, bandwidth){
  print("---------------------------")
  print(model)
  tas_model <- tas_cmip5[tas_cmip5$model == model, ]
  tas_model_factual <- tas_model[tas_model$experiment == "historical" | tas_model$experiment == "rcp85", ]
  print(head(tas_model_factual))
  if(model == "HadCRUT"){
    tas_model_counterfactual <- tas_model_factual[tas_model_factual$year <= 1900,]
  } else {
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
  }
  print(head(tas_model_counterfactual))
  # tas_model_counterfactual <- tas_model[tas_model$experiment == "historical" & tas_model$year <= 1900, ]
  p12 <- estim_p12_ns(x = tas_model_counterfactual$tas, t = tas_model_factual$year, z = tas_model_factual$tas, kernel = kernel, h = bandwidth)
}, tas_cmip5 = tas_cmip5, kernel = krnl, bandwidth = bandwidth)

lp12boot_ggplot <- mapply(
  function(p12, model){
    data.frame(
      model = model,
      year = p12$t_unique,
      p12 = p12$p12_hat,
      ci_q05 = p12$p12_hat - qnorm(0.05) * p12$sigma_p12_hat, 
      ci_q95 = p12$p12_hat + qnorm(0.05) * p12$sigma_p12_hat 
    )
  },  p12 = lp12, model = lmodels, SIMPLIFY = FALSE
) %>% do.call(rbind, .)
p <- ggplot(lp12boot_ggplot) +
  geom_hline(yintercept = 1/2) +
  geom_line(aes(y = p12, x = year, colour = model), lwd = 1.2) +
  geom_ribbon(aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = model), alpha = 0.3) +
  facet_wrap( ~ model, ncol = 5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +  ylim(min(lp12boot_ggplot$ci_q95), 1) + 
  ggtitle("q1(t), factual = historical + rcp85, counterfactual = historicalNat") + ylab("q1")
plot(p)
ggsave("p12_ggplot.pdf", dpi = "retina", width = 20, height = 15, units = "cm")


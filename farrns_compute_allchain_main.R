library(ncdf4)
library(magrittr)
library(ggplot2)
devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("~/2_Code/Naveau/farFAR/farr/inprogress/tests_nsfar_algo.R")
source("check_H0_algo.R")
source("farrns_compute_allchain_algo.R")

tas_cmip5 <- readRDS(file = "tas_cmip5.rds")
nbtimestep <- aggregate(tas ~ model + experiment, data = tas_cmip5, FUN = length)

lmodels <- unique(tas_cmip5$model) %>% as.character()
lmodels <- subset(nbtimestep, experiment == "rcp85", select = model) %>% unlist() %>% as.character()

tas_cmip5 <- tas_cmip5[tas_cmip5$model %in% lmodels, ]

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


tas_cmip5 <- rbind(tas_cmip5, tas_hadcrut)
lmodels <- unique(tas_cmip5$model) %>% as.character()


kernel <- kernel_epanechnikov
tas_hadcrut_counterfactual <- subset(tas_hadcrut, year <= 1900)
theta <- estim_theta.nswexp(x = tas_hadcrut_counterfactual$tas, t = tas_hadcrut$year, z = tas_hadcrut$tas, kernel = kernel, h = NULL)
bandwidth <- theta$h


p12_multimodel <- compute_allchain(lmodels, tas_cmip5, kernel, bandwidth)


lp12boot_ggplot <- mapply(
  function(p12, model){
    data.frame(
      model = model,
      year = p12$tpred,
      p12 = p12$p12_hat,
      ci_q05 = p12$p12_hat - qnorm(0.95) * p12$sigma_p12_hat, 
      ci_q95 = p12$p12_hat + qnorm(0.95) * p12$sigma_p12_hat 
    )
  },  p12 = p12_multimodel$lp12, model = lmodels, SIMPLIFY = FALSE
) %>% do.call(rbind, .)
p12_allinone_ggplot <- rbind(
  lp12boot_ggplot,
  cbind(
    model = "Multi-Model",
    p12_multimodel$p12_multimodel[, 1:2],
    ci_q05 = p12_multimodel$p12_multimodel[, 2] - qnorm(0.95) * p12_multimodel$p12_multimodel[, 3],
    ci_q95 = p12_multimodel$p12_multimodel[, 2] + qnorm(0.95) * p12_multimodel$p12_multimodel[, 3]
  )
)

p <- ggplot(p12_allinone_ggplot) +
  geom_hline(yintercept = 1/2) +
  geom_line(aes(y = p12, x = year, colour = model), lwd = 1.2) +
  geom_ribbon(aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = model), alpha = 0.3) +
  facet_wrap( ~ model, ncol = 5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +  ylim(min(p12_allinone_ggplot$ci_q05, na.rm =TRUE), 1) + 
  ggtitle("q1(t), factual = historical + rcp85, counterfactual = historicalNat") + ylab("q1")
plot(p)

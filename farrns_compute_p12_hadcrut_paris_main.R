library(ncdf4)
library(magrittr)
library(ggplot2)
library(maps)
library(reshape2)
library(fields)
source("tests_nsfar_algo.R")
source("farrns_compute_allchain_AD_algo.R")
source("check_H0_simpler_algo.R")

# Load Data -----
hadcrut_nc <- "tas_hadcrut_augustavg.nc"
nc <- nc_open(file = hadcrut_nc)
lon <- ncvar_get(nc, "longitude") %>% as.numeric()
lat <- ncvar_get(nc, "latitude") %>% as.numeric()
year <- ncvar_get(nc, "time") %/% 10000 %>% as.numeric()
tas <- ncvar_get(nc, "temperature_anomaly")
tas <- t(matrix(tas, ncol = length(year)))
nc_close(nc)
tas_hadcrut <- data.frame(
  "institute" = "HadCRUT",
  "model" = "HadCRUT",
  "experiment" = "historical",
  "run" = "obs",
  "year" = year
)
tas_hadcrut <- cbind(tas_hadcrut, tas)
tas_hadcrut_paris <-
  tas_hadcrut[
    c("institute", "model", "experiment", "run", "year", as.character(iparis))
  ]
tas_hadcrut_paris_counterfactual <- subset(tas_hadcrut_paris, year <= 1900)


# p12 Analysis
tpred <- min(tas_hadcrut_paris$year):max(tas_hadcrut_paris$year)
kernel <- kernel_epanechnikov
bandwidth <- 32
p12_hadcrut_paris <- estim_p12_ns(
  x = tas_hadcrut_paris_counterfactual[, as.character(iparis)],
  t = tas_hadcrut_paris$year,
  z = tas_hadcrut_paris[, as.character(iparis)],
  tpred = tpred, kernel = kernel, h = bandwidth
)

# Reformat for plots
hadcrut_paris_plot <- with(
  p12_hadcrut_paris,
  data.frame(
    year = tpred,
    p12 = p12_hat,
    ci_q05 = p12_hat - qnorm(0.95) * sigma_p12_hat,
    ci_q95 = p12_hat + qnorm(0.95) * sigma_p12_hat
  )
)


pdf(
  paste0("tas_hadcrut_augustavg_paris.pdf"),
  width = 20 / 2.54, height = 20 / 2.54
)
p <- ggplot(hadcrut_paris_plot) +
  geom_hline(yintercept = 1 / 2, lwd = 0.1) +
  geom_line(aes(y = p12, x = year), colour = "grey", lwd = 0.5) +
  geom_ribbon(
    aes(ymin = ci_q05, ymax = ci_q95, x = year),
    colour = "grey", alpha = 0.3
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(
    paste0(
      varname_inplot,
      ", q(t), HadCRUT, factual = 1850-2018, counterfactual = 1850-1900"
    )
  ) +
  ylab("q(t)")
plot(p)
dev.off()


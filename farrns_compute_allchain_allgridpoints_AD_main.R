library(ncdf4)
library(magrittr)
library(ggplot2)
library(maps)
library(reshape2)
library(fields)
library(viridis)
library(metR)
# devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("tests_nsfar_algo.R")
source("farrns_compute_allchain_AD_algo.R")
source("check_H0_simpler_algo.R")

# Analysis parameters -----
hadcrut_nc <- "tas_hadcrut_augustavg.nc"

# cmip5_rds <- "tas_cmip5_augustavg.rds"
# cmip5_prefix <- paste0(
#   sub(pattern = "(.*)\\..*$", replacement = "\\1", cmip5_rds),
#   "_onerun_experts"
# )
# cmip6_rds <- "tas_cmip6_augustavg.rds"
# p12outputs_rds <- "p12_augustavg_cmip56_experts_onerun.rds"
# varname_inplot <- "tas, august mean"

cmip5_rds <- "pr_cmip5_yearmax.rds"
cmip5_prefix <- paste0(
  sub(pattern = "(.*)\\..*$", replacement = "\\1", cmip5_rds),
  "_onerun_experts"
)
cmip6_rds <- "pr_cmip6_yearmax.rds"
p12outputs_rds <- "p12_pryearmax_cmip56_experts_onerun.rds"
varname_inplot <- "pr, yearmax"


# Data loading and preparation ----
worldmap <- maps::map("world", plot = FALSE)

nc <- nc_open(file = hadcrut_nc)
lon <- ncvar_get(nc, "longitude") %>% as.numeric()
lat <- ncvar_get(nc, "latitude") %>% as.numeric()
year <- ncvar_get(nc, "time") %/% 10000 %>% as.numeric()
tas <- ncvar_get(nc, "temperature_anomaly")
image(
  lon[lon >= -10 & lon <= 30],
  lat[lat >= 35 & lat <= 70],
  tas[lon >= -10 & lon <= 30, lat >= 35 & lat <= 70, 15]
)
lines(worldmap)
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

tas_cmip5 <- readRDS(file = cmip5_rds)
tas_cmip6 <- readRDS(file = cmip6_rds)
tas_cmip6 <- subset(tas_cmip6, !(model == "CanESM5(CM6)" & run == "r1i1p1f1"))
# aggregate(run ~ model + institute, data = tas_cmip6, FUN = function(x) paste(unique(x), collapse = ", "))
tas_cmip5 <- rbind(tas_cmip5, tas_cmip6)
# levels(tas_cmip5$model) <- toupper(levels(tas_cmip5$model))
# levels(tas_cmip5$institute) <- toupper(levels(tas_cmip5$institute))
nbtimestep <- aggregate(
  tas_cmip5$"1" ~ model + experiment,
  data = tas_cmip5, FUN = length
)
names(nbtimestep)[3] <- "nbtimestep"
lmodels <- names(table(nbtimestep$model))[table(nbtimestep$model) == 3]
nbtimestep_historicalNat <- subset(nbtimestep, experiment == "historicalNat" &
  nbtimestep >= 100)
firstyear <- aggregate(year ~ model + experiment, data = tas_cmip5, FUN = min)
firstyear <- aggregate(
  year ~ model,
  data = subset(firstyear, experiment %in% c("historical", "historicalNat")),
  FUN = function(x) all(x == 1850)
)
firstyear <- subset(firstyear, year == TRUE)

lmodels <- names(table(nbtimestep$model))[table(nbtimestep$model) == 3]
lmodels <- intersect(lmodels, nbtimestep_historicalNat$model)
lmodels <- intersect(lmodels, firstyear$model)
tas_cmip5 <- tas_cmip5[tas_cmip5$model %in% lmodels, ]
run_table <- aggregate(
  run ~ model,
  data = tas_cmip5,
  FUN = function(x) paste(unique(x)[1])
)
tas_cmip5 <- mapply(
  function(model, run) {
    tas_cmip5[tas_cmip5$model == model & tas_cmip5$run == run, ]
  },
  model = run_table$model, run = run_table$run,
  SIMPLIFY = FALSE
) %>%
  do.call(rbind, .)

lonlat_df <- expand.grid(lon = lon, lat = lat)
iworld <- seq.int(nrow(lonlat_df))
iparis <- with(lonlat_df, which(lon == 2.5 & lat == 47.5))

lmodels <- unique(tas_cmip5[, c("model", "institute")])$model %>%
  as.character()
linstitutes <- unique(tas_cmip5[, c("model", "institute")])$institute %>%
  as.character()

# Hadcrut analysis ----
tas_hadcrut_paris <-
  tas_hadcrut[
    c("institute", "model", "experiment", "run", "year", as.character(iparis))
  ]
tas_hadcrut_paris_counterfactual <- subset(tas_hadcrut_paris, year <= 1900)
tpred <- min(tas_hadcrut_paris$year):max(tas_hadcrut_paris$year)
kernel <- kernel_epanechnikov
bandwidth <- 32
p12_hadcrut_paris <- estim_p12_ns(
  x = tas_hadcrut_paris_counterfactual[, as.character(iparis)],
  t = tas_hadcrut_paris$year,
  z = tas_hadcrut_paris[, as.character(iparis)],
  tpred = tpred, kernel = kernel, h = bandwidth
)

hadcrut_paris_plot <- with(
  p12_hadcrut_paris,
  data.frame(
    year = tpred,
    p12 = p12_hat,
    ci_q05 = p12_hat - qnorm(0.95) * sigma_p12_hat,
    ci_q95 = p12_hat + qnorm(0.95) * sigma_p12_hat
  )
)


# Multimodel analysis ----

lp12_pergridpoint <- lapply(iworld, function(igridpoint) {
  cat(".")
  tas_cmip5 <-
    tas_cmip5[
      ,
      c("institute", "model", "experiment", "run", "year", igridpoint)
    ]
  names(tas_cmip5)[6] <- "tas"
  ina <- is.na(tas_cmip5$tas)
  tas_cmip5 <- tas_cmip5[!ina, ]
  # head(tas_cmip5) %>% print()
  p12_multimodel <- compute_allchain(
    lmodels, linstitutes, tas_cmip5, kernel, bandwidth
  )
})
saveRDS(lp12_pergridpoint, file = p12outputs_rds)
lp12_pergridpoint <- readRDS(p12outputs_rds)

reformat_lp12 <- function(lp12_pergridpoint, lmodels, linstitutes) {
  p12 <- mapply(function(p12_onegridpoint, igridpoint) {
    p12_multimodel_df <- data.frame(
      igridpoint = igridpoint,
      institute = "multimodel",
      model = "multimodel",
      year = p12_onegridpoint$p12_multimodel$year,
      p12 = p12_onegridpoint$p12_multimodel$p12,
      sig = p12_onegridpoint$p12_multimodel$sig
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
    rbind(lp12_model_df, p12_multimodel_df)
  },
  p12_onegridpoint = lp12_pergridpoint,
  igridpoint = iworld,
  SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
}
lp12_reformated <- reformat_lp12(
  lp12_pergridpoint = lp12_pergridpoint,
  lmodels = lmodels, linstitutes = linstitutes
)

p12_paris <- subset(lp12_reformated, igridpoint == iparis)
p12_paris$ci_q05 <- p12_paris$p12 - qnorm(0.95) * p12_paris$sig
p12_paris$ci_q95 <- p12_paris$p12 + qnorm(0.95) * p12_paris$sig

# Plots ----
pdf(paste0(cmip5_prefix, "_paris.pdf"), width = 20 / 2.54, height = 20 / 2.54)
# png(paste0(cmip5_prefix,"_paris.png"), units = "in", width = 20/2.54, height = 20/2.54, res = 150)
p <- ggplot(subset(p12_paris, model != "multimodel")) +
  geom_hline(yintercept = 1 / 2, lwd = 0.1) +
  geom_line(aes(y = p12, x = year, colour = institute), lwd = 0.5) +
  geom_ribbon(
    aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = institute),
    alpha = 0.3
  ) +
  facet_wrap(~ institute + model, ncol = 5) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(
    paste0(
      varname_inplot,
      ", q1(t), factual = historical + rcp85, counterfactual = historicalNat"
    )
  ) +
  ylab("q1")
plot(p)
dev.off()

pdf(
  paste0(cmip5_prefix, "_multimodel_paris.pdf"),
  width = 20 / 2.54, height = 20 / 2.54
)
p <- ggplot(subset(p12_paris, model == "multimodel")) +
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
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(
    paste0(
      varname_inplot,
      ", q1(t), factual = historical + rcp85, counterfactual = historicalNat"
    )
  ) +
  ylab("q1")
plot(p)
dev.off()
pdf(
  paste0(cmip5_prefix, "_hadcrut_paris.pdf"),
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
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(
    paste0(
      varname_inplot,
      ", q1(t), HadCRUT, factual = 1850-2018, counterfactual = 1850-1900"
    )
  ) +
  ylab("q1")
plot(p)
dev.off()
pdf(
  paste0(cmip5_prefix, "_multimodel_paris.pdf"),
  width = 20 / 2.54, height = 20 / 2.54
)
p <- ggplot(subset(p12_paris, model == "multimodel")) +
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
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(
    paste0(
      varname_inplot,
      ", q1(t), factual = historical + rcp85, counterfactual = historicalNat"
    )
  ) +
  ylab("q1")
plot(p)
dev.off()
pdf(
  paste0(cmip5_prefix, "_hadcrut_vs_multimodel_paris.pdf"),
  width = 25 / 2.54, height = 20 / 2.54
)
p <- ggplot(
  rbind(
    cbind(model = "HadCRUT", hadcrut_paris_plot),
    subset(
      p12_paris,
      model == "multimodel",
      select = c(model, year, p12, ci_q05, ci_q95)
    )
  )
) +
  geom_hline(yintercept = 1 / 2, lwd = 0.1) +
  geom_line(aes(y = p12, x = year, colour = model, group = model), lwd = 2) +
  geom_ribbon(
    aes(ymin = ci_q05, ymax = ci_q95, x = year, group = model, fill = model),
    alpha = 0.3
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  # xlim (1850, 2100) + ggtitle("p12 (t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(paste0(varname_inplot, ", q1(t), HadCRUT vs multimodel, factual = 1850-2018, counterfactual = 1850-1900")) +
  ylab("q1")
plot(p)
dev.off()

pdf(paste0(cmip5_prefix, "_maps_multimodel.pdf"), width = 20 / 2.54, height = 20 / 2.54)
p12_4years <- subset(
  lp12_reformated,
  year %in% round(seq.int(from = 1850, to = 2100, length.out = 9)) &
    model == "multimodel"
)
p12_4years$lon <- lonlat_df$lon[as.numeric(paste(p12_4years$igridpoint))]


p12_4years$lat <- lonlat_df$lat[as.numeric(paste(p12_4years$igridpoint))]
WorldData <- map_data("world")
colpal <- rev(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))
p <- ggplot(p12_4years) +
  #   geom_raster(aes(x = lon, y = lat, fill = p12)) +
  geom_contour_fill(aes(x = lon, y = lat, z = p12), breaks = c(-0.01, seq(0.1, 0.9, by = 0.1), 1.01)) +
  geom_map(
    data = WorldData, map = WorldData,
    aes(map_id = region), fill = NA, col = "black"
  ) +
  facet_wrap(~year, ncol = 3) +
  #   coord_cartesian(xlim=c(-10,30), ylim=c(35, 70)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradientn(colours = colpal, limits = c(-0.01, 1.01), breaks = c(0, seq(0.1, 0.9, by = 0.1), 1)) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historic al[1850-1900]")
  ggtitle(paste0(varname_inplot, ", q1(t) multimodel, factual = historical + rcp85, counterfactual = historicalNat"))
plot(p)

dev.off()
lp12_paris <- lp12_pergridpoint[[which(iworld == iparis)]]
# CvM_paris_filtered <- keep_onemodel_perinstitute(
#   lp12_paris$CvM_df,
#   unique(tas_cmip5[, c("institute", "model")])
# )
# CvM_paris_filtered$weight <-
#   CvM_paris_filtered$weight / sum(CvM_paris_filtered$weight)
# CvM_paris_filtered <- within(
#   CvM_paris_filtered,
#   model <- factor(model, levels = model[order(weight, decreasing = TRUE)])
# )
CvM_paris_filtered <- lp12_paris$CvM_df
CvM_paris_filtered <- melt(
  CvM_paris_filtered,
  id.vars = c("institute", "model")
)
pdf(paste0(cmip5_prefix, "_CvM.pdf"), width = 20 / 2.54, height = 20 / 2.54)
p <- ggplot(
  data = CvM_paris_filtered,
  aes(x = model, y = value, fill = grepl("(CM6)", model))
) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(variable), scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(
    values = c("darkgrey", "red"),
    labels = c("CMIP5", "CMIP6"),
    name = "CMIP"
  )
plot(p)
dev.off()

# tas_hadcrut_counterfactual_paris <- subset(tas_cmip5, model == "HadCRUT" & year <= 1900)[[paste(iparis)]]
pdf(paste0(cmip5_prefix, "_hist_checkH0.pdf"), width = 20 / 2.54, height = 20
/ 2.54)
p <- hist_checkH0(
  lUhat = lp12_paris$lUhat,
  lmodels = lmodels,
  linstitutes = linstitutes
) +
  theme(
    legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)
  )
plot(p)
dev.off()

pdf(
  paste0(cmip5_prefix, "_qqplot_checkH0.pdf"),
  width = 20 / 2.54, height = 25 / 2.54
)
p <- qqplot_checkH0(
  lUhat = lp12_paris$lUhat,
  lmodels = lmodels,
  linstitutes = linstitutes
) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
plot(p)
dev.off()
# Question de philippes sur les dates des max de Uhat
# mapply(function(indice, model){subset(tas_cmip5, model == model & year <=
# 1900)[indice, "year"]}, indice = sapply(lp12_paris$lUhat[lmodels !=
# "HadCRUT"], which.max), model = lmodels[lmodels != "HadCRUT"])
gc(TRUE)

CvM_allgridpoints <- mapply(
  function(x, igridpoint) {
    cbind(x$CvM, gridpoint = as.numeric(igridpoint))
  },
  x = lp12_pergridpoint, igridpoint = iworld,
  SIMPLIFY = FALSE
) %>% do.call(rbind, .)
CvM_allgridpoints$gridpoint <- apply(
  lonlat_df[CvM_allgridpoints$gridpoint, ],
  1,
  function(x) paste0("[", x[1], ";", x[2], "]")
)
CvM_allgridpoints <- within(
  CvM_allgridpoints,
  gridpoint <- factor(gridpoint, levels = unique(gridpoint)),
  model <- factor(model, levels = lmodels[order(linstitutes)])
)
pdf(paste0(cmip5_prefix, "_CvM_allgripoints.pdf"), width = 20 / 2.54, height = 25 / 2.54)
ylim <- c(0, with(CvM_allgridpoints, max(CvM[is.finite(CvM)])))
p <- ggplot(subset(CvM_allgridpoints)) +
  geom_raster(aes(x = model, y = gridpoint, fill = CvM)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_fill_gradientn(colours = rev(magma(5)), limits = ylim) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  ggtitle(paste0("Checking H0: CvM distances for each location and each model"))
plot(p)
dev.off()
CvM_matrix <- acast(
  data = CvM_allgridpoints,
  formula = model ~ gridpoint,
  value.var = "CvM"
)
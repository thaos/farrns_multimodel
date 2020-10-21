library(ncdf4)
library(magrittr)
library(ggplot2)
library(maps)
library(reshape2)
library(fields)
library(metR)
library(viridis)
# devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("tests_nsfar_algo.R")
source("farrns_compute_allchain_AD_algo.R")
source("check_H0_simpler_algo.R")

# ---------------------------------
# Parameters -----
# ---------------------------------

config <- new.env()
source("config_pryearmax.R",  local = config) 
cmip5_prefix <- config$cmip5_prefix
p12outputs_rds <- config$p12outputs_rds 
varname_inplot <- config$varname_inplot
varunit <- config$varunit
if(!dir.exists(cmip5_prefix)) dir.create(cmip5_prefix)

lmodels <- readRDS(
  paste0(cmip5_prefix, "_lmodels.rds")
)
linstitutes <- readRDS(
  paste0(cmip5_prefix, "_linstitutes.rds")
)

hadcrut_nc <- "tas_hadcrut_augustavg.nc"
nc <- nc_open(file = hadcrut_nc)
lon <- ncvar_get(nc, "longitude") %>% as.numeric()
lat <- ncvar_get(nc, "latitude") %>% as.numeric()
year <- ncvar_get(nc, "time") %/% 10000 %>% as.numeric()
nc_close(nc)

lonlat_df <- expand.grid(lon = lon, lat = lat)
iworld <- seq.int(nrow(lonlat_df))
iparis <- with(lonlat_df, which(lon == 2.5 & lat == 47.5))
ineighbours <- with(
  lonlat_df,
  which(
    lon == 2.5 & lat == (47.5 + 5) |
    lon == 2.5 & lat == (47.5 - 5) |
    lon == (2.5 + 5) & lat == 47.5 |
    lon == (2.5 - 5) & lat == 47.5
  ) 
)

# ---------------------------------
# Loadings Results ---
# ---------------------------------

linstitutes  <- readRDS(
  paste0(cmip5_prefix, "_linstitutes.rds")
)
linstitutes  <- readRDS(
  paste0(cmip5_prefix, "_linstitutes.rds")
)
lp12_pergridpoint <- readRDS(p12outputs_rds)
lp12_reformated <- reformat_lp12(
  lp12_pergridpoint = lp12_pergridpoint,
  lmodels = lmodels, linstitutes = linstitutes
)


# ---------------------------------
# Plots ----
# ---------------------------------


# ---------------------------------
# q(t) Paris plots

p12_paris <- subset(lp12_reformated, igridpoint == iparis)
p12_paris$ci_q05 <- p12_paris$p12 - qnorm(0.95) * p12_paris$sig
p12_paris$ci_q95 <- p12_paris$p12 + qnorm(0.95) * p12_paris$sig

pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_paris.pdf")),
  width = 20 / 2.54, height = 20 / 2.54
)
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
      ", q(t), factual = historical + rcp85, counterfactual = historicalNat"
    )
  ) +
  ylab("q(t)")
plot(p)
dev.off()

pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_multimodel_paris.pdf")),
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
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(
    paste0(
      varname_inplot,
      ", q(t), factual = historical + rcp85, counterfactual = historicalNat"
    )
  ) +
  ylab("q(t)")
plot(p)
dev.off()


# ---------------------------------
## q(t) worldmaps

p12_4years <- subset(
  lp12_reformated,
  year %in% c(1850, 1900, 1940, 1970, 2000, 2020, 2030, 2050, 2100) &
    model == "multimodel"
)
p12_4years$lon <- lonlat_df$lon[as.numeric(paste(p12_4years$igridpoint))]


p12_4years$lat <- lonlat_df$lat[as.numeric(paste(p12_4years$igridpoint))]

WorldData <- map_data("world")
colpal <- rev(
  c(
    "#9e0142", "#d53e4f", "#f46d43",
    "#fdae61", "#fee08b", "#ffffbf",
    "#e6f598", "#abdda4", "#66c2a5",
    "#3288bd", "#5e4fa2"
  )
)
pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_maps_multimodel.pdf")),
  width = 20 / 2.54, height = 20 / 2.54
)
p <- ggplot(p12_4years) +
  geom_contour_fill(aes(x = lon, y = lat, z = p12), breaks = c(-0.01, seq(0.1, 0.9, by = 0.1), 1.01)) +
  geom_map(
    data = WorldData, map = WorldData,
    aes(map_id = region), fill = NA, col = "black"
  ) +
  facet_wrap(~year, ncol = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradientn(name = "q(t)", colours = colpal, limits = c(-0.01, 1.01), breaks = c(0, seq(0.1, 0.9, by = 0.1), 1)) +
  ggtitle(paste0(varname_inplot, ", q(t) multimodel, factual = historical + rcp85, counterfactual = historicalNat"))
plot(p)
dev.off()


# ---------------------------------
## Checking Assumption
lp12_paris <- lp12_pergridpoint[[iparis]]
pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_qqplot_checkH0.pdf")),
  width = 20 / 2.54, height = 25 / 2.54
)
p <- qqplot_checkH0(
  lXmZm = lp12_paris$lXmZm,
  lmodels = lmodels,
  linstitutes = linstitutes
) +
  xlab(paste0("Xm (", varunit, ")")) +
  ylab(paste0("Zm (", varunit, ")")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
plot(p)
dev.off()


# ---------------------------------
## Expert Aggregation

### Weights in Paris

lp12_paris <- lp12_pergridpoint[[which(iworld == iparis)]]
# )
weights_paris <- data.frame(
  model = factor(lmodels, levels = lmodels[order(linstitutes)]),
  weight = lp12_paris$lweights
)
pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_weights.pdf")), 
  width = 20 / 2.54, height = 20 / 2.54
)
p <- ggplot(
  data = weights_paris,
  aes(x = model, y = weight, fill = grepl("(CM6)", model))
) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(
    values = c("darkgrey", "red"),
    labels = c("CMIP5", "CMIP6"),
    name = "CMIP"
  )
plot(p)
dev.off()

### Weights in Paris neighbours
weights_neighbours <- mapply(
  function(x, igridpoint) {
    data.frame(model = lmodels, weight = x$lweights, gridpoint = as.numeric(igridpoint))
  },
  x = lp12_pergridpoint[ineighbours], igridpoint = ineighbours,
  SIMPLIFY = FALSE
) %>% do.call(rbind, .)
weights_neighbours$gridpoint <- apply(
  lonlat_df[weights_neighbours$gridpoint, ],
  1,
  function(x) paste0("[", x[1], ";", x[2], "]")
)
weights_neighbours <- within(
  weights_neighbours,
  gridpoint <- factor(gridpoint, levels = unique(gridpoint)),
  model <- factor(lmodel, levels = lmodels[order(linstitutes)])
)
ylim <- c(0, with(weights_neighbours, max(weight[is.finite(weight)])))
pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_weights_neighbours.pdf")),
  width = 20 / 2.54, height = 25 / 2.54
)
p <- ggplot(
  data  = weights_neighbours,
  aes(x = model, y = weight, fill = grepl("(CM6)", model))
  ) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(
    values = c("darkgrey", "red"),
    labels = c("CMIP5", "CMIP6"),
    name = "CMIP"
  ) +
  facet_wrap(~ gridpoint, ncol = 2) +
  ggtitle(paste0("Expert aggregation: weights distances for gridpoints around Paris"))
plot(p)
dev.off()

### Weights all gridpoints
weights_allgridpoints <- mapply(
  function(x, igridpoint) {
    data.frame(model = lmodels, weight = x$lweights, gridpoint = as.numeric(igridpoint))
  },
  x = lp12_pergridpoint, igridpoint = iworld,
  SIMPLIFY = FALSE
) %>% do.call(rbind, .)
weights_allgridpoints$gridpoint <- apply(
  lonlat_df[weights_allgridpoints$gridpoint, ],
  1,
  function(x) paste0("[", x[1], ";", x[2], "]")
)
weights_allgridpoints <- within(
  weights_allgridpoints,
  gridpoint <- factor(gridpoint, levels = unique(gridpoint)),
  model <- factor(lmodel, levels = lmodels[order(linstitutes)])
)
ylim <- c(0, with(weights_allgridpoints, max(weight[is.finite(weight)])))

pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_weights_allgripoints.pdf")),
  width = 20 / 2.54, height = 25 / 2.54
)
p <- ggplot(weights_allgridpoints) +
  geom_raster(aes(x = model, y = gridpoint, fill = weight)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_fill_gradientn(colours = rev(magma(5)), limits = ylim) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  ggtitle(paste0("Expert aggregation: weights distances for each location and each model"))
plot(p)
dev.off()




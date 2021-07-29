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
obs_lonlat <- config$obs_lonlat
p12obs_rds <- config$p12obs_rds
city <- config$city

if(!dir.exists(cmip5_prefix)) dir.create(cmip5_prefix)


hadcrut_nc <- "tas_hadcrut_augustavg.nc"
nc <- nc_open(file = hadcrut_nc)
lon <- ncvar_get(nc, "longitude") %>% as.numeric()
lat <- ncvar_get(nc, "latitude") %>% as.numeric()
year <- ncvar_get(nc, "time") %/% 10000 %>% as.numeric()
nc_close(nc)

lonlat_df <- expand.grid(lon = lon, lat = lat)
iworld <- seq.int(nrow(lonlat_df))
iparis <- with(lonlat_df, which(lon == obs_lonlat[1] & lat == obs_lonlat[2]))
ineighbours <- with(
  lonlat_df,
  which(
    lon == obs_lonlat[1] & lat == (obs_lonlat[2] + 5) |
    lon == obs_lonlat[1] & lat == (obs_lonlat[2] - 5) |
    lon == (obs_lonlat[1] + 5) & lat == obs_lonlat[2] |
    lon == (obs_lonlat[1] - 5) & lat == obs_lonlat[2]
  ) 
)

# ---------------------------------
# Loadings Results ---
# ---------------------------------

lmodels <- readRDS(
  paste0(cmip5_prefix, "_lmodels.rds")
)
linstitutes <- readRDS(
  paste0(cmip5_prefix, "_linstitutes.rds")
)
lp12_pergridpoint <- readRDS(p12outputs_rds)
p12_obs <- readRDS(p12obs_rds)
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
weight_paris_df <- data.frame(
  institute = linstitutes,
  model = lmodels,
  weight = lp12_pergridpoint[[iparis]]$lweights_kl
)

GmZ_paris <- mapply(
    function(institute, model, lp12){
    data.frame(institute = institute, model = model, year = lp12$tz, GmZ =lp12$GmZ, stringsAsFactors = FALSE)},
    institute = linstitutes, model = lmodels, lp12 = lp12_pergridpoint[[iparis]]$lp12,
    SIMPLIFY = FALSE
)
GmZ_paris <- do.call(rbind, GmZ_paris)

pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_", city, ".pdf")),
  width = 20 / 2.54, height = 20 / 2.54
)
p <- ggplot(GmZ_paris) +
  geom_point(
    aes(x = year,
        y = GmZ
    ),
    size = 0.1
  ) +
  geom_hline(yintercept = 1 / 2, lwd = 0.1) +
  geom_ribbon(
    data = subset(p12_paris, !(model %in% c("multimodel_kl" , "multimodel_oracle", "multimodel_best"))),
    aes(ymin = ci_q05, ymax = ci_q95, x = year),
    fill = "orange",
    alpha = 0.8
  ) +
  geom_line(
    data = subset(p12_paris, !(model %in% c("multimodel_kl" , "multimodel_oracle", "multimodel_best"))),
    aes(y = p12, x = year),
    colour = "red",
    lwd = 1
  ) +
  #geom_point(
  #  data = weight_paris_df,
  #  aes(x = 1860,
  #      y = 0.9,
  #      fill = weight, size = weight
  #  ),
  #  shape = 23
  #) +
  geom_label(
      data = weight_paris_df,
      aes(x = 1875,
          y = 0.9,
          label = sprintf("%1.2f", weight)
      ),
      size = 3
  ) + 
  scale_fill_gradientn(colours = rev(magma(5)), limits = c(0, 1)) +
  facet_wrap(~ institute + model, ncol = ceiling(sqrt(length(lmodels)))) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  guides(size = "none") +
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

p12_obs_df <- data.frame(
  igridpoint = iparis, 
  institute = "obs",
  model = paste(city, "obs"),
  year = p12_obs$tpred,
  p12 = p12_obs$p12_hat,
  sig = p12_obs$sigma_p12_hat,
  ci_q05 = p12_obs$p12_hat - qnorm(0.95) * p12_obs$sigma_p12_hat,
  ci_q95 = p12_obs$p12_hat + qnorm(0.95) * p12_obs$sigma_p12_hat
)

pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_multimodel_", city, ".pdf")),
  width = 20 / 2.54, height = 20 / 2.54
)
mm_df <- subset(
    p12_paris,
    model == "multimodel_kl" | model == "multimodel_best"
)
mm_df$combination <- factor(mm_df$model, level = c("multimodel_best", "multimodel_kl"))
levels(mm_df$combination)[levels(mm_df$combination)=="multimodel_best"] <- "GCM with largest Anderson-Darling p-value"
levels(mm_df$combination)[levels(mm_df$combination)=="multimodel_kl"] <- "Weighted convex combination"
p <- ggplot(mm_df)+
  geom_hline(yintercept = 1 / 2, lwd = 0.1) +
  geom_ribbon(
    aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = combination),
    alpha = 0.3
  ) +
  geom_line(aes(y = p12, x = year,  colour = combination), lwd = 2) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  guides(fill=guide_legend(title=NULL), colour=guide_legend(title=NULL)) +
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(
    paste0(
      varname_inplot,
      ", q(t), factual = hist + rcp85, counterfactual = histNat"
    )
  ) +
  ylab("q(t)")
plot(p)
dev.off()



pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_obs_", city, ".pdf")),
  width = 20 / 2.54, height = 20 / 2.54
)
p <- ggplot(
  rbind(
    p12_obs_df,
    subset(
      p12_paris,
      model == "multimodel_kl" 
    )
  )) +
  geom_hline(yintercept = 1 / 2, lwd = 0.1) +
  geom_line(aes(y = p12, x = year, colour = model), lwd = 0.5) +
  geom_ribbon(
    aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = model),
    alpha = 0.3
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  xlim(1850, 2100) +
  ylim(min(p12_paris$ci_q05), 1) +
  ggtitle(
    paste0(
      varname_inplot,
      ", q(t), factual = hist + rcp85, counterfactual = histNat"
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
    model == "multimodel_kl"
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
  width = 20 / 2.54, height = 15 / 2.54
)
p <- ggplot(p12_4years) +
  geom_contour_fill(aes(x = lon, y = lat, z = p12), breaks = c(-0.01, seq(0.1, 0.9, by = 0.1), 1.01)) +
  geom_map(
    data = WorldData, map = WorldData,
    aes(map_id = region), fill = NA, col = "black"
  ) +
  facet_wrap(~year, ncol = 3) +
  theme(
      # axis.text.x = element_text(angle = 90, hjust = 1),
      legend.key.height = unit(0.7, "inch"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white')
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(name = "q(t)", colours = colpal, limits = c(-0.01, 1.01), breaks = c(0, seq(0.1, 0.9, by = 0.1), 1)) +
  coord_fixed() +
  ggtitle(paste0(varname_inplot, ", q(t) multimodel"))
plot(p)
dev.off()

p12_reformated_klmm <- subset(lp12_reformated, model == "multimodel_kl" & year <= 2100)
emergence_df <- by(p12_reformated_klmm, p12_reformated_klmm$igridpoint, function(df){
  ci_q05 <- df$p12 - qnorm(0.95) * df$sig
  ci_q95 <- df$p12 + qnorm(0.95) * df$sig
  year <- df$year
  emergence <- sign <- NA
  icond <- ci_q05 > 0.5
  i <- length(icond)
  while(i > 0 & icond[i]){
    emergence <- year[i]
    i <- i - 1
  }
  if(!is.na(emergence)){
    sign <- "+"
    return(cbind(df[1,], emergence = emergence, sign = sign))
  }
  icond <- ci_q95 < 0.5
  i <- length(icond)
  while(i > 0 & icond[i]){
    emergence <- year[i]
    i <- i - 1
  }
  if(!is.na(emergence)){
    sign <- "-"
    return(cbind(df[1,], emergence = emergence, sign = sign))
  }
  return(cbind(df[1,], emergence = emergence, sign = sign))
  }
) %>% do.call(rbind, .)

emergence_df$lon <- lonlat_df$lon[emergence_df$igridpoint]
emergence_df$lat <- lonlat_df$lat[emergence_df$igridpoint]


pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_emergence_multimodel.pdf")),
  width = 20 / 2.54, height = 12 / 2.54
)
zlim <- c(1900, 2100)
p <- ggplot(subset(emergence_df, sign == "+"))+
  geom_raster(aes(x = lon, y = lat, fill = emergence)) +
  geom_map(
    data = WorldData, map = WorldData,
    aes(map_id = region), fill = NA, col = "black"
  ) +
  theme(
      legend.key.height = unit(0.5, "inch"),
      # axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white')
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(name = "year", colours = colpal,  breaks = seq(zlim[1], zlim[2], by = 20), limits = zlim) +
  coord_fixed() +
  scale_shape_manual(values=c(24, 25))+
  ggtitle(paste0(varname_inplot, ", multimodel time of emergence"))
plot(p)
dev.off()

# ---------------------------------
## Checking Assumption
lp12_paris <- lp12_pergridpoint[[iparis]]
pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_qqplot_checkH0_", city, ".pdf")),
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
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
plot(p)
dev.off()


# ---------------------------------
## Expert Aggregation

### Weights in Paris

pdf(
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_weights.pdf_", city, ".pdf")), 
  width = 20 / 2.54, height = 20 / 2.54
)
p <- ggplot(
  data = weight_paris_df,
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
    data.frame(model = lmodels, weight = x$lweights_kl, gridpoint = as.numeric(igridpoint))
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
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_weights_neighbours_", city, ".pdf")),
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
  ggtitle(paste0("Aggregation: KL weights for gridpoints around ", city))
plot(p)
dev.off()

### Weights all gridpoints
weights_allgridpoints <- mapply(
  function(x, igridpoint) {
    data.frame(model = lmodels, weight = x$lweights_kl, gridpoint = as.numeric(igridpoint))
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

write.table(
  aggregate(weight ~ model, weights_allgridpoints,  function(w) mean(w == 0)),
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_weights_summary_mean.txt"))
)
write.table(
  aggregate(weight ~ model, weights_allgridpoints,  function(w) mean(w)),
  file.path(cmip5_prefix, paste0(cmip5_prefix, "_weights_summary_zerofreq.txt"))
)

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
  ggtitle(paste0("Expert aggregation: weights for each location and each model"))
plot(p)
dev.off()




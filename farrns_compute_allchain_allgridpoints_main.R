library(ncdf4)
library(magrittr)
library(ggplot2)
library(maps)
library(reshape2)
source("tests_nsfar_algo.R")
source("check_H0_algo.R")
source("farrns_compute_allchain_algo.R")

# Parameters --------------------------------------------------------------
cmip5_rds <- "tas_cmip5_augustavg.rds"
cmip5_prefix <- sub(pattern = "(.*)\\..*$", replacement = "\\1", cmip5_rds)
cmip6_rds <- "tas_cmip6_augustavg.rds"
hadcrut_nc <- "tas_hadcrut_augustavg.nc"
p12outputs_rds <- "p12_augustavg_cmip56.rds"
varname_inplot <- "tas, august mean"


worldmap <- maps::map("world", plot = FALSE)

tas_df <- readRDS("tas_df.rds")
# tas_df a data.frame with columns institute, model, experiement, run, year, 
# and the average temperature in august for different grid points. 
# For those, the name of the colomn corresponds to the row number of the data.frame
lonlat_df <- readRDS("lonlat_df.rds")
# lonlat_df which contains the coordinates of each HadCRUT grid point for the whole world.
# lonlat_df a data.frame with colunms lon and lat which gives the coordinates of the grid points.
# the ith row of the data.frame gives the coordinates of the grid point named "i".

# iparis, the gridpoint name for the grid point closest to Paris  
iparis <- with(lonlat_df, which(lon == 2.5 & lat == 47.5))

# list of grid point indices for gridpoints without NA values in the HadCRUT dataset.
inotna_gridpoint <- names(tas_df)[-(1:5)]
lmodels <- unique(tas_df[, c("model", "institute")])$model %>% as.character()
linstitutes <- unique(tas_df[, c("model", "institute")])$institute%>% as.character()

kernel <- kernel_epanechnikov
bandwidth <- 32

# compute p12(t) for all models and for the multimodel
lp12_pergridpoint <- lapply(inotna_gridpoint, function(igridpoint){
  tas_df <- tas_df[, c("institute", "model", "experiment", "run", "year", igridpoint)]
  names(tas_df)[6] <- "tas"
  ina <- is.na(tas_df$tas)
  tas_df <- tas_df[!ina, ]
  head(tas_df) %>% print()
  p12_multimodel <- compute_allchain(lmodels, tas_df, kernel, bandwidth)
})
saveRDS(lp12_pergridpoint, file = p12outputs_rds)

# reformat results into a data.frame
reformat_lp12 <- function(lp12_pergridpoint, lmodels, linstitutes, inotna_gridpoint){
  p12 <- mapply(function(p12_onegridpoint, igridpoint){
    p12_multimodel_df <- data.frame(
      igridpoint =  igridpoint,
      institute = "multimodel",
      model = "multimodel",
      year = p12_onegridpoint$p12_multimodel$year,
      p12 = p12_onegridpoint$p12_multimodel$p12,
      sig = p12_onegridpoint$p12_multimodel$sig
      )
    lp12_model_df <- mapply(function(p12_df, model, institute){
      p12_model_df <- data.frame(
        igridpoint =  igridpoint,
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
igridpoint = inotna_gridpoint,
SIMPLIFY = FALSE
) %>% do.call(rbind, .)
}
lp12_reformated <- reformat_lp12(
  lp12_pergridpoint = lp12_pergridpoint,
  lmodels = lmodels, linstitutes = linstitutes,
  inotna_gridpoint = inotna_gridpoint
  )

p12_paris <- subset(lp12_reformated, igridpoint == iparis)
p12_paris$ci_q05 = p12_paris$p12 - qnorm(0.95) * p12_paris$sig 
p12_paris$ci_q95 = p12_paris$p12 + qnorm(0.95) * p12_paris$sig 

pdf(paste0(cmip5_prefix,"_paris.pdf"), width = 20/2.54, height = 20/2.54)
p <- ggplot(p12_paris) +
  geom_hline(yintercept = 1/2, lwd = 0.1) +
  geom_line(aes(y = p12, x = year, colour = institute), lwd = 0.5) +
  geom_ribbon(aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = institute), alpha = 0.3) +
  facet_wrap( ~ institute + model, ncol = 5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +  ylim(min(p12_paris$ci_q05), 1) + 
  ggtitle(paste0(varname_inplot, ", q1(t), factual = historical + rcp85, counterfactual = historicalNat")) +
  ylab("q1")
plot(p)
dev.off()

# pdf(paste0(cmip5_prefix,"_maps.pdf"), width = 20/2.54, height = 15/2.54)
# for(year_toplot in seq(1850, 2100, by = 10)){
#   p12_oneyear <- subset(lp12_reformated, year == year_toplot)
#   p12_oneyear$lon <- lonlat_df$lon[as.numeric(paste(p12_oneyear$igridpoint))]
#   p12_oneyear$lat <- lonlat_df$lat[as.numeric(paste(p12_oneyear$igridpoint))]
#   
#   WorldData <- map_data('world')
#   colpal <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
#   p <- ggplot(p12_oneyear) +
#     geom_raster(aes(x = lon, y = lat, fill = p12)) +
#     geom_map(
#       data = WorldData, map = WorldData,
#       aes(map_id = region), fill = NA, col = "black"
#     ) +
#     facet_wrap( ~ model, ncol = 5) +
#     coord_cartesian(xlim=c(-10,30), ylim=c(35, 70)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     scale_fill_gradientn(colours = colpal, limits = c(0, 1), values = seq(0, 1, length.out = 11), breaks = seq(0, 1, length.out = 11)) +
#     # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
#     ggtitle(paste0(varname_inplot,", q1(", year_toplot, "), factual = historical + rcp85, counterfactual = historicalNat"))
#   plot(p)
# }
# dev.off()
pdf(paste0(cmip5_prefix,"_maps_multimodel.pdf"), width = 20/2.54, height = 20/2.54)
p12_4years <- subset(lp12_reformated, year %in% round(seq.int(from = 1850, to = 2100, length.out = 9)) & model == "multimodel")
p12_4years$lon <- lonlat_df$lon[as.numeric(paste(p12_4years$igridpoint))]
p12_4years$lat <- lonlat_df$lat[as.numeric(paste(p12_4years$igridpoint))]
WorldData <- map_data('world')
colpal <- rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'))
p <- ggplot(p12_4years) +
  geom_raster(aes(x = lon, y = lat, fill = p12)) +
  geom_map(
    data = WorldData, map = WorldData,
    aes(map_id = region), fill = NA, col = "black"
    ) +
  facet_wrap( ~ year, ncol = 3) +
  coord_cartesian(xlim=c(-10,30), ylim=c(35, 70)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradientn(colours = colpal, limits = c(0, 1), values = seq(0, 1, length.out = 11), breaks = seq(0, 1, length.out = 11)) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  ggtitle(paste0(varname_inplot,", q1(t) multimodel, factual = historical + rcp85, counterfactual = historicalNat"))
plot(p)
dev.off()
lp12_paris <- lp12_pergridpoint[[which(inotna_gridpoint == iparis)]]
CvM_paris_filtered <- keep_onemodel_perinstitute(lp12_paris$CvM_df, unique(tas_df[, c("institute", "model")]))
CvM_paris_filtered <- subset(CvM_paris_filtered, model != "HadCRUT")
CvM_paris_filtered$weight   <- CvM_paris_filtered$weight /  sum(CvM_paris_filtered$weight)
CvM_paris_filtered <- within(
  CvM_paris_filtered, 
  model <- factor(model,  levels=model[order(weight, decreasing = TRUE)])
  )
CvM_paris_filtered <- melt(CvM_paris_filtered, id.vars = c("institute", "model"))
pdf(paste0(cmip5_prefix,"_CvM.pdf"), width = 20/2.54, height = 20/2.54)
p <- ggplot(data = CvM_paris_filtered, aes(x = model, y = value, fill = grepl("(CM6)", model))) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(variable),  scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("darkgrey", "red"), labels = c("CMIP5", "CMIP6"), name= "CMIP")
plot(p)
dev.off()  

tas_hadcrut_counterfactual_paris <- subset(tas_df, model == "HadCRUT" & year <= 1900)[[paste(iparis)]]
pdf(paste0(cmip5_prefix,"_hist_checkH0.pdf"), width = 20/2.54, height = 20/2.54)
p <- hist_checkH0(
  X = tas_hadcrut_counterfactual_paris,
  lUhat = lp12_paris$lUhat[lmodels != "HadCRUT"],
  lmodels = lmodels[lmodels != "HadCRUT"],
  linstitutes = linstitutes[linstitutes != "HadCRUT"]
  ) +
theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
  plot(p)
dev.off() 

pdf(paste0(cmip5_prefix,"_qqplot_checkH0.pdf"), width = 20/2.54, height = 24/2.54)
p <- qqplot_checkH0(
  X = tas_hadcrut_counterfactual_paris,
  lUhat = lp12_paris$lUhat[lmodels != "HadCRUT"],
  lmodels = lmodels[lmodels != "HadCRUT"],
  linstitutes = linstitutes[linstitutes != "HadCRUT"]
  ) 
plot(p)
dev.off() 
# Question de philippes sur les dates des max de Uhat
# mapply(function(indice, model){subset(tas_df, model == model & year <= 1900)[indice, "year"]}, indice = sapply(lp12_paris$lUhat[lmodels != "HadCRUT"], which.max), model = lmodels[lmodels != "HadCRUT"])

library(ncdf4)
library(magrittr)
library(ggplot2)
library(maps)
library(reshape2)
library(fields)
library(viridis)
library(metR)
library(ggplot2)
library(gganimate)

# devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("tests_nsfar_algo.R")
source("farrns_compute_allchain_AD_algo.R")
source("check_H0_simpler_algo.R")

# Parameters --------------------------------------------------------------
cmip5_rds <- "tas_cmip5_augustavg.rds"
cmip5_prefix <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", cmip5_rds), "_onerun_CvMsimple")
cmip6_rds <- "tas_cmip6_augustavg.rds"
hadcrut_nc <- "tas_hadcrut_augustavg.nc"
p12outputs_rds <- "p12_augustavg_cmip56_CvMsimple_onerun.rds"
varname_inplot <- "tas, august mean"

# cmip5_rds <- "pr_cmip5_yearmax.rds"
# cmip5_prefix <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", cmip5_rds), "_onerun_CvMsimple")
# cmip6_rds <- "pr_cmip6_yearmax.rds"
# p12outputs_rds <- "p12_pryearmax_cmip56_CvMsimple_onerun.rds"
# varname_inplot <- "pr, yearmax"

worldmap <- maps::map("world", plot = FALSE)

nc = nc_open(file = hadcrut_nc)
lon = ncvar_get(nc, "longitude") %>% as.numeric()
lat = ncvar_get(nc, "latitude") %>% as.numeric()
year = ncvar_get(nc, "time") %/% 10000 %>% as.numeric()
tas = ncvar_get(nc, "temperature_anomaly")
image(lon[lon >= - 10 & lon <= 30], lat[lat >= 35 & lat <= 70], tas[lon >= - 10 & lon <= 30, lat >= 35 & lat <= 70, 15])
lines(worldmap)
tas = t(matrix(tas, ncol = length(year)))
nc_close(nc)
tas_hadcrut <- data.frame("institute" = "HadCRUT", "model" = "HadCRUT", "experiment" = "historical", "run" = "obs", "year" = year)
tas_hadcrut = cbind(tas_hadcrut, tas)

tas_cmip5 <- readRDS(file = cmip5_rds)
tas_cmip6 <- readRDS(file = cmip6_rds)
tas_cmip6 <- subset(tas_cmip6, !(model == "CanESM5(CM6)" & run ==  "r1i1p1f1"))
tas_cmip5 <- rbind(tas_cmip5, tas_cmip6)
nbtimestep <- aggregate(tas_cmip5$"1" ~ model + experiment, data = tas_cmip5, FUN = length)
names(nbtimestep)[3] <- "nbtimestep"
lmodels <- names(table(nbtimestep$model))[table(nbtimestep$model) == 3]
nbtimestep_historicalNat <- subset(nbtimestep, experiment == "historicalNat" &
                                     nbtimestep >= 100)  
lmodels <- intersect(lmodels, nbtimestep_historicalNat$model)
tas_cmip5 <- tas_cmip5[tas_cmip5$model %in% lmodels, ]
run_table <- aggregate(run ~ model, data = tas_cmip5, FUN = function(x) paste(unique(x)[1]))
tas_cmip5 <- mapply(
  function(model, run){
    tas_cmip5[tas_cmip5$model == model & tas_cmip5$run == run, ]
  },
  model = run_table$model, run = run_table$run,
  SIMPLIFY = FALSE
) %>% 
  do.call(rbind, .)

lonlat_df <- expand.grid(lon = lon, lat = lat)
iworld <- seq.int(nrow(lonlat_df))
iparis <- with(lonlat_df, which(lon == 2.5 & lat == 47.5))

lmodels <- unique(tas_cmip5[, c("model", "institute")])$model %>% as.character()
linstitutes <- unique(tas_cmip5[, c("model", "institute")])$institute %>% as.character()
rm(tas_cmip5, tas_cmip6, tas_hadcrut)

reformat_lp12 <- function(lp12_pergridpoint, lmodels, linstitutes){
  p12 <- mapply(function(p12_onegridpoint, igridpoint){
    p12_multimodel_df <- data.frame(
      igridpoint =  igridpoint,
      institute = "multimodel",
      model = "multimodel",
      year = p12_onegridpoint$p12_multimodel$year,
      p12 = p12_onegridpoint$p12_multimodel$p12,
      sig = p12_onegridpoint$p12_multimodel$sig
    )
    p12_multimodel_df
  }, 
  p12_onegridpoint = lp12_pergridpoint, 
  igridpoint = iworld,
  SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
}
lp12_pergridpoint <- readRDS(p12outputs_rds)
lp12_reformated <- reformat_lp12(
  lp12_pergridpoint = lp12_pergridpoint,
  lmodels = lmodels, linstitutes = linstitutes
)

p12_allyears <- subset(lp12_reformated, year %in% seq.int(from = 1850, to = 2100) & model == "multimodel")
p12_allyears$lon <- lonlat_df$lon[as.numeric(paste(p12_allyears$igridpoint))]
p12_allyears$lat <- lonlat_df$lat[as.numeric(paste(p12_allyears$igridpoint))]
rm(lp12_reformated, lp12_pergridpoint)
WorldData <- map_data('world')
colpal <- rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'))

resfactor = 3
for(year in 1850:2100){
  print(year)
  p12_year = p12_allyears[p12_allyears$year == year, ]
  filename = sprintf("anim_img/%s_%05d.png", cmip5_prefix, year)
  png(filename,  height = resfactor * 480, width = resfactor * 720, res = resfactor * 72)
  p <- ggplot(p12_year) +
    #   geom_raster(aes(x = lon, y = lat, fill = p12)) +
    geom_contour_fill(aes(x = lon, y = lat, z = p12), breaks = c(-0.01, seq(0.1, 0.9, by = 0.05), 1.01)) +
    geom_map(
      data = WorldData, map = WorldData,
      aes(map_id = region), fill = NA, col = "black"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradientn(colours = colpal, limits = c(-0.01, 1.01), breaks = c(0, seq(0.1, 0.9, by = 0.1), 1)) +
    ggtitle(paste0(varname_inplot, ", ", year,', p1(t) multimodel, factual = historical + rcp85, counterfactual = historicalNat'))
    # Here comes the gganimate specific bits
    # transition_time(year) +
    # ease_aes('linear')
  plot(p)
  dev.off()
}
  

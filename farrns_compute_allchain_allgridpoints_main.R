library(ncdf4)
library(magrittr)
library(ggplot2)
library(maps)
library(reshape2)
devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("~/2_Code/Naveau/farFAR/farr/inprogress/tests_nsfar_algo.R")
source("check_H0_algo.R")
source("farrns_compute_allchain_algo.R")

for(iconfig in 1:1){
  
  if(iconfig == 1) {
    # Parameters --------------------------------------------------------------
    cmip5_rds <- "tas_cmip5_augustavg.rds"
    cmip5_prefix <- sub(pattern = "(.*)\\..*$", replacement = "\\1", cmip5_rds)
    cmip6_rds <- "tas_cmip6_augustavg.rds"
    hadcrut_nc <- "tas_hadcrut_augustavg.nc"
    p12outputs_rds <- "p12_augustavg_cmip56.rds"
    varname_inplot <- "tas, august mean"
  }
  
  
  if(iconfig == 2) {
    cmip5_rds <- "tas_cmip5_yearmax.rds"
    hadcrut_nc <- "tas_hadcrut_yearmax.nc"
    p12outputs_rds <- "p12_yearmax.rds"
    varname_inplot <- "tas, year max"
  }
  
  if(iconfig == 3) {
    cmip5_rds <- "tas_cmip5_yearmin.rds"
    hadcrut_nc <- "tas_hadcrut_yearmin.nc"
    p12outputs_rds <- "p12_yearmin.rds"
    varname_inplot <- "tas, year min"
  } 
  # End parameters
  
  
  worldmap <- maps::map("world", plot = FALSE)
  
  tas_cmip5 <- readRDS(file = cmip5_rds)
  tas_cmip6 <- readRDS(file = cmip6_rds)
  tas_cmip6 <- subset(tas_cmip6, !(model == "CanESM5(CM6)" & run ==  "r1i1p1f1"))
  # aggregate(run ~ model + institute, data = tas_cmip6, FUN = function(x) paste(unique(x), collapse = ", "))
  tas_cmip5 <- rbind(tas_cmip5, tas_cmip6)
  # levels(tas_cmip5$model) <- toupper(levels(tas_cmip5$model)) 
  # levels(tas_cmip5$institute) <- toupper(levels(tas_cmip5$institute)) 
  nbtimestep <- aggregate(tas_cmip5$"1" ~ model + experiment, data = tas_cmip5, FUN = length)
  lmodels <- names(table(nbtimestep$model))[table(nbtimestep$model) == 3]
  
  tas_cmip5 <- tas_cmip5[tas_cmip5$model %in% lmodels, ]
  
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
  
  lonlat_df <- expand.grid(lon = lon, lat = lat)
  ieurope <- with(lonlat_df, which(lon >= - 10 & lon <= 30 & lat >= 35 & lat <= 70))
  iparis <- with(lonlat_df, which(lon == 2.5 & lat == 47.5))
  tas_hadcrut <- tas_hadcrut[, c("institute", "model", "experiment", "run", "year", paste(ieurope))]
  tas_cmip5 <- tas_cmip5[, c("institute", "model", "experiment", "run", "year", paste(ieurope))]
  
  ina_gridpoint <- apply(tas_hadcrut, 2, function(x) sum(is.na(x)) > 10)
  tas_hadcrut <- tas_hadcrut[, !ina_gridpoint]
  tas_cmip5 <- tas_cmip5[, !ina_gridpoint]
  inotna_gridpoint <- names(tas_hadcrut)[-(1:5)]
  
  tas_cmip5 <- rbind(tas_cmip5, tas_hadcrut)
  lmodels <- unique(tas_cmip5$model) %>% as.character()
  if(cmip5_rds == "tas_cmip5_yearmax.rds"){
    tas_cmip5 <- aggregate(data = tas_cmip5, . ~  year + run + experiment + model + institute, FUN = max)
  }
  if(cmip5_rds == "tas_cmip5_yearmin.rds"){
    tas_cmip5 <- aggregate(data = tas_cmip5, . ~  year + run + experiment + model + institute, FUN = min)
  }
  
  kernel <- kernel_epanechnikov
  bandwidth <- 32
  
  lp12_pergridpoint <- lapply(inotna_gridpoint, function(igridpoint){
    tas_cmip5 <- rbind(
      tas_cmip5[, c("institute", "model", "experiment", "run", "year", igridpoint)],
      tas_hadcrut[, c("institute", "model", "experiment", "run", "year", igridpoint)]
    )
    names(tas_cmip5)[6] <- "tas"
    ina <- is.na(tas_cmip5$tas)
    tas_cmip5 <- tas_cmip5[!ina, ]
    head(tas_cmip5) %>% print()
    p12_multimodel <- compute_allchain(lmodels, tas_cmip5, kernel, bandwidth)
  })
  saveRDS(lp12_pergridpoint, file = p12outputs_rds)
  
  reformat_lp12 <- function(lp12_pergridpoint, lmodels, inotna_gridpoint){
    p12 <- mapply(function(p12_onegridpoint, igridpoint){
      p12_multimodel_df <- data.frame(
        igridpoint =  igridpoint,
        model = "multimodel",
        year = p12_onegridpoint$p12_multimodel$year,
        p12 = p12_onegridpoint$p12_multimodel$p12,
        sig = p12_onegridpoint$p12_multimodel$sig
      )
      lp12_model_df <- mapply(function(p12_df, model){
        p12_model_df <- data.frame(
          igridpoint =  igridpoint,
          model = model,
          year = p12_df$tpred,
          p12 = p12_df$p12_hat,
          sig = p12_df$sigma_p12_hat
        )
      },
      p12_df = p12_onegridpoint$lp12, model = lmodels,
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
    lmodels = lmodels,
    inotna_gridpoint = inotna_gridpoint
  )
  
  p12_paris <- subset(lp12_reformated, igridpoint == iparis)
  p12_paris$ci_q05 = p12_paris$p12 - qnorm(0.95) * p12_paris$sig 
  p12_paris$ci_q95 = p12_paris$p12 + qnorm(0.95) * p12_paris$sig 
  
  pdf(paste0(cmip5_prefix,"_paris.pdf"), width = 20/2.54, height = 15/2.54)
  p <- ggplot(p12_paris) +
    geom_hline(yintercept = 1/2, lwd = 0.1) +
    geom_line(aes(y = p12, x = year, colour = model), lwd = 0.5) +
    geom_ribbon(aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = model), alpha = 0.3) +
    facet_wrap( ~ model, ncol = 5) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
    xlim(1850, 2100) +  ylim(min(p12_paris$ci_q05), 1) + 
    ggtitle(paste0(varname_inplot, ", q1(t), factual = historical + rcp85, counterfactual = historicalNat")) +
    ylab("q1")
  plot(p)
  dev.off()
  
  pdf(paste0(cmip5_prefix,"_maps.pdf"), width = 20/2.54, height = 15/2.54)
  for(year_toplot in seq(1850, 2100, by = 10)){
    p12_oneyear <- subset(lp12_reformated, year == year_toplot)
    p12_oneyear$lon <- lonlat_df$lon[as.numeric(paste(p12_oneyear$igridpoint))]
    p12_oneyear$lat <- lonlat_df$lat[as.numeric(paste(p12_oneyear$igridpoint))]
    
    WorldData <- map_data('world')
    colpal <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
    p <- ggplot(p12_oneyear) +
      geom_raster(aes(x = lon, y = lat, fill = p12)) +
      geom_map(
        data = WorldData, map = WorldData,
        aes(map_id = region), fill = NA, col = "black"
      ) +
      facet_wrap( ~ model, ncol = 5) +
      coord_cartesian(xlim=c(-10,30), ylim=c(35, 70)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradientn(colours = colpal, limits = c(0, 1), values = seq(0, 1, length.out = 11), breaks = seq(0, 1, length.out = 11)) +
      # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
      ggtitle(paste0(varname_inplot,", q1(", year_toplot, "), factual = historical + rcp85, counterfactual = historicalNat"))
    plot(p)
  }
  dev.off()
  pdf(paste0(cmip5_prefix,"_maps_multimodel.pdf"), width = 20/2.54, height = 15/2.54)
  p12_4years <- subset(lp12_reformated, year %in% round(seq.int(from = 1850, to = 2100, length.out = 9)) & model == "multimodel")
  p12_4years$lon <- lonlat_df$lon[as.numeric(paste(p12_4years$igridpoint))]
  p12_4years$lat <- lonlat_df$lat[as.numeric(paste(p12_4years$igridpoint))]
  WorldData <- map_data('world')
  colpal <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
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
  CvM_paris_filtered <- keep_onemodel_perinstitute(lp12_paris$CvM_df, unique(tas_cmip5[, c("institute", "model")]))
  CvM_paris_filtered <- subset(CvM_paris_filtered, model != "HadCRUT")
  CvM_paris_filtered$weight   <- CvM_paris_filtered$weight /  sum(CvM_paris_filtered$weight)
  CvM_paris_filtered <- within(
    CvM_paris_filtered, 
    model <- factor(model,  levels=model[order(weight, decreasing = TRUE)])
  )
  CvM_paris_filtered <- melt(CvM_paris_filtered, id.vars = c("institute", "model"))
  pdf(paste0(cmip5_prefix,"_CvM.pdf"), width = 20/2.54, height = 15/2.54)
  p <- ggplot(data = CvM_paris_filtered, aes(x = model, y = value)) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(variable),  scales = "free") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
  plot(p)
  dev.off()  
  
  tas_hadcrut_counterfactual_paris <- subset(tas_cmip5, model == "HadCRUT" & year <= 1900)[[paste(iparis)]]
  pdf(paste0(cmip5_prefix,"_hist_checkH0.pdf"), width = 20/2.54, height = 15/2.54)
   p <- hist_checkH0(
    X = tas_hadcrut_counterfactual_paris,
    lUhat = lp12_paris$lUhat[lmodels != "HadCRUT"],
    lmodels = lmodels[lmodels != "HadCRUT"]
  ) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
  plot(p)
  dev.off() 
  
  pdf(paste0(cmip5_prefix,"_qqplot_checkH0.pdf"), width = 20/2.54, height = 15/2.54)
  p <- qqplot_checkH0(
    X = tas_hadcrut_counterfactual_paris,
    lUhat = lp12_paris$lUhat[lmodels != "HadCRUT"],
    lmodels = lmodels[lmodels != "HadCRUT"]
  ) 
  plot(p)
  dev.off()  
  gc(TRUE)
}

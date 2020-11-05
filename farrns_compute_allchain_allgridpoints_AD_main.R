library(ncdf4)
library(magrittr)
library(reshape2)
source("tests_nsfar_algo.R")
source("farrns_compute_allchain_AD_algo.R")
source("check_H0_simpler_algo.R")

# ---------------------------------
# Parameters -----
# ---------------------------------

config <- new.env()
source("config_tasaugustavg.R",  local = config) 
cmip5_rds <- config$cmip5_rds
cmip6_rds <- config$cmip6_rds
cmip5_prefix <- config$cmip5_prefix
p12outputs_rds <- config$p12outputs_rds 
unit_scaling <- config$unit_scaling
obs_rds <- config$obs_rds
obs_lonlat <- config$obs_lonlat
obs_varname <- config$obs_varname
p12obs_rds <- config$p12obs_rds

# ---------------------------------
# Data loading and preparation ----
# ---------------------------------

tas_cmip5 <- readRDS(file = cmip5_rds)
tas_cmip6 <- readRDS(file = cmip6_rds)
tas_cmip6 <- subset(tas_cmip6, !(model == "CanESM5(CM6)" & run == "r1i1p1f1"))
tas_cmip5 <- rbind(tas_cmip5, tas_cmip6)


# ---------------------------------
# First Filter: Only complete models and runs

## Keep only models that have the 3 experiments
nbtimestep <- aggregate(
  tas_cmip5$"1" ~ model + experiment,
  data = tas_cmip5, FUN = length
)
names(nbtimestep)[3] <- "nbtimestep"
lmodels <- names(table(nbtimestep$model))[table(nbtimestep$model) == 3]

## Check that first year is historical and historicalNat run is 1850
firstyear <- aggregate(year ~ model + experiment, data = tas_cmip5, FUN = min)
firstyear <- aggregate(
  year ~ model,
  data = subset(firstyear, experiment %in% c("historical", "historicalNat")),
  FUN = function(x) all(x == 1850)
)
firstyear <- subset(firstyear, year == TRUE)
lmodels <- intersect(lmodels, firstyear$model)

## Keep only data for the corresponding models
tas_cmip5 <- tas_cmip5[tas_cmip5$model %in% lmodels, ]


# ---------------------------------
# Second Filter: keep only the first run

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

# Update remaining list of models and institutes
lmodels <- unique(tas_cmip5[, c("model", "institute")])$model %>%
  as.character()
linstitutes <- unique(tas_cmip5[, c("model", "institute")])$institute %>%
  as.character()
saveRDS(
  lmodels,
  paste0(cmip5_prefix, "_lmodels.rds")
)
saveRDS(
  linstitutes,
  paste0(cmip5_prefix, "_linstitutes.rds")
)


# ---------------------------------
# Rescale the variables
tas_cmip5[, -(1:5)] <- unit_scaling(tas_cmip5[, -(1:5)])

# ---------------------------------
# Hadcrut lat lon
hadcrut_nc <- "tas_hadcrut_augustavg.nc"
nc <- nc_open(file = hadcrut_nc)
lon <- ncvar_get(nc, "longitude") %>% as.numeric()
lat <- ncvar_get(nc, "latitude") %>% as.numeric()
year <- ncvar_get(nc, "time") %/% 10000 %>% as.numeric()
nc_close(nc)

lonlat_df <- expand.grid(lon = lon, lat = lat)
iworld <- seq.int(nrow(lonlat_df))
iparis <- with(lonlat_df, which(lon == obs_lonlat[1] & lat == obs_lonlat[2]))

# ---------------------------------
# Multimodel analysis ----
# ---------------------------------
compute_allbandwidth <- function(lmodels, linstitutes, tas_cmip5, kernel) {
  lbandwidth <- lapply(lmodels, function(model, tas_cmip5, kernel) {
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" | tas_model$experiment == "rcp85", ]
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
    tpred <- min(tas_model$year):max(tas_model$year)
    p12 <- cross_validate_h(x = tas_model_counterfactual$tas, t = tas_model_factual$year, z = tas_model_factual$tas, tpred = tpred, kernel = kernel)
  }, tas_cmip5 = tas_cmip5, kernel = kernel)
  return(unlist(lbandwidth))
}
tas_paris <-
	tas_cmip5[
		  ,
		  c("institute", "model", "experiment", "run", "year", iparis)
		  ]
names(tas_paris)[6] <- "tas"
ina <- is.na(tas_paris$tas)
tas_paris <- tas_paris[!ina, ]

kernel <- kernel_epanechnikov
allbandwidths <- compute_allbandwidth(
				   lmodels, linstitutes, tas_paris, kernel
				   )
print(str(allbandwidths))
bandwidth <- median(allbandwidths)
print(bandwidth)

obs_df <- readRDS(obs_rds)
p12_obs <- estim_p12_ns(x = obs_df[, obs_varname][obs_df[, "year"] >= 1850 & obs_df[, "year"] <= 1900], t = obs_df[, "year"], z = obs_df[, obs_varname], tpred = obs_df[, "year"], kernel = kernel, h = bandwidth)
saveRDS(p12_obs, file = p12obs_rds)

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
  p12_multimodel <- compute_allchain(
    lmodels, linstitutes, tas_cmip5, kernel, bandwidth
  )
})
saveRDS(lp12_pergridpoint, file = p12outputs_rds)
lp12_pergridpoint <- readRDS(p12outputs_rds)




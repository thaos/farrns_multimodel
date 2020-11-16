cmip5_rds <- "pr_cmip5_yearmax.rds"
cmip6_rds <- "pr_cmip6_yearmax.rds"
# obs_rds <- "pr_hohenpeissenberg.rds"
obs_rds <- "pr_ghcdn.rds"
obs_lonlat <- c(-2.5, 52.5)
# obs_lonlat <- c(-27.5, -62.5) #min(toe)
# obs_lonlat <- c(12.5, 47.5)
obs_varname <- "pr"
p12outputs_rds <- "p12_pryearmax_cmip56_experts_onerun.rds"
p12obs_rds <- "p12_pryearmax_obs.rds"
cmip5_prefix <-  "pr_cmip5_yearmax_onerun_experts"
varname_inplot <- "pr, yearmax"
unit_scaling <- function(x) x * 3600 * 24
varunit <- "mm/day"
city <- "oxford"
# city <- "-27.5E_-62.5N"
# city <- "hohenpeissenberg"



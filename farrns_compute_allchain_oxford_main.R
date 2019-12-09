library(ncdf4)
library(magrittr)
library(ggplot2)
library(maps)
library(reshape2)
# devtools::load_all("~/2_Code/Naveau/farFAR/farr")
source("~/2_Code/Naveau/farFAR/farr/inprogress/tests_nsfar_algo.R")
source("check_H0_algo.R")
source("farrns_compute_allchain_algo.R")

# Parameters --------------------------------------------------------------
cmip5_rds <- "pr_cmip5_yearmax.rds"
cmip5_prefix <- sub(pattern = "(.*)\\..*$", replacement = "\\1", cmip5_rds)
cmip6_rds <- "pr_cmip6_yearmax.rds"
ghcdn_rds <- "pr_ghcdn.rds"
p12outputs_rds <- "p12_sonmax_cmip56.rds"
varname_inplot <- "pr, SON max"
  
worldmap <- maps::map("world", plot = FALSE)

pr_cmip5 <- readRDS(file = cmip5_rds)
pr_cmip6 <- readRDS(file = cmip6_rds)
pr_cmip6 <- subset(pr_cmip6, !(model == "CanESM5(CM6)" & run ==  "r1i1p1f1"))
# aggregate(run ~ model + institute, data = pr_cmip6, FUN = function(x) paste(unique(x), collapse = ", "))
pr_cmip5 <- rbind(pr_cmip5, pr_cmip6)
pr_cmip5$pr <- pr_cmip5$pr * 24 * 3600
# levels(pr_cmip5$model) <- toupper(levels(pr_cmip5$model)) 
# levels(pr_cmip5$institute) <- toupper(levels(pr_cmip5$institute)) 
nbtimestep <- aggregate(pr_cmip5$pr ~ model + experiment, data = pr_cmip5, FUN = length)
lmodels <- names(table(nbtimestep$model))[table(nbtimestep$model) == 3]

pr_cmip5 <- pr_cmip5[pr_cmip5$model %in% lmodels, ]
run_check <- aggregate(run ~ experiment + model + institute, data = pr_cmip5, FUN = function(x) unique(c(paste(x))))
run_check <- dcast(
  data = run_check, 
  model + institute ~ experiment
)
for(i in 1:nrow(run_check)){
  print("-----------")
  run_model <- run_check[i, ]
  allruns <- union(union(run_model$historicalNat[[1]], run_model$historical[[1]]), run_model$rcp85[[1]])
  run_tokeep <- intersect(intersect(run_model$historicalNat[[1]], run_model$historical[[1]]), run_model$rcp85[[1]])
  run_torm <- setdiff(allruns, run_tokeep)
  print(run_model$model[[1]] )
  print(run_torm)
  if(length(run_torm) > 0){
    pr_cmip5 <- pr_cmip5[!(pr_cmip5$model == paste(run_model$model[[1]]) & pr_cmip5$run %in% run_torm), ]
  }
  print(dim(pr_cmip5))
}

pr_ghcnd <- data.frame("institute" = "HadCRUT", "model" = "HadCRUT", "experiment" = "historical", "run" = "obs")
pr_ghcnd = cbind(pr_ghcnd, readRDS(ghcdn_rds))

pr_cmip5 <- rbind(pr_cmip5, pr_ghcnd)
pr_cmip5 <- pr_cmip5[!is.na(pr_cmip5$pr), ]
pr_cmip5 <- pr_cmip5[!(pr_cmip5$model == "BNU-ESM"), ]
lmodels <- unique(pr_cmip5[, c("model", "institute")])$model %>% as.character()
linstitutes <- unique(pr_cmip5[, c("model", "institute")])$institute%>% as.character()


kernel <- kernel_epanechnikov
# theta <- estim_theta.nswexp(x = pr_ghcnd$pr[pr_ghcnd$year <= 1900], t = pr_ghcnd$year, z = pr_ghcnd$pr, kernel = kernel, h = NULL)
# bandwidth <- theta$h
bandwidth <- 31

pr_cmip5 <- rbind(pr_cmip5, pr_ghcnd)
pr_cmip5$tas <- pr_cmip5$pr
p12_allchain <- compute_allchain(lmodels, pr_cmip5, kernel, bandwidth)

reformat_lp12 <- function(p12_allchain, lmodels, linstitutes){
  p12_multimodel_df <- cbind(
    institute = "multimodel",
    model = "multimodel",
    p12_allchain$p12_multimodel
  )
  lp12_model_df <- mapply(function(p12_df, model, institute){
    p12_model_df <- data.frame(
      institute = institute,
      model = model,
      year = p12_df$tpred,
      p12 = p12_df$p12_hat,
      sig = p12_df$sigma_p12_hat
    )
  },  p12_df = p12_allchain$lp12, model = lmodels, institute = linstitutes,
  SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  rbind(lp12_model_df, p12_multimodel_df)
}
lp12_reformated <- reformat_lp12(
  p12_allchain = p12_allchain,
  lmodels = lmodels, linstitutes = linstitutes
)
# 
lp12_reformated$ci_q05 = lp12_reformated$p12 - qnorm(0.95) * lp12_reformated$sig 
lp12_reformated$ci_q95 = lp12_reformated$p12 + qnorm(0.95) * lp12_reformated$sig 
levels(lp12_reformated$institute)[levels(lp12_reformated$institute) ==  "HadCRUT"] <- "GHCDN"
levels(lp12_reformated$model)[levels(lp12_reformated$model) ==  "HadCRUT"] <- "GHCDN"
# 
pdf(paste0(cmip5_prefix,"_oxford.pdf"), width = 20/2.54, height = 20/2.54)
p <- ggplot(lp12_reformated) +
  geom_hline(yintercept = 1/2, lwd = 0.1) +
  geom_line(aes(y = p12, x = year, colour = institute), lwd = 0.5) +
  geom_ribbon(aes(ymin = ci_q05, ymax = ci_q95, x = year, fill = institute), alpha = 0.3) +
  facet_wrap( ~ institute + model, ncol = 5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  # xlim(1850, 2100) + ggtitle("p12(t), counterfactual = historical[1850-1900]")
  xlim(1850, 2100) +  ylim(min(lp12_reformated$ci_q05), 1) +
  ggtitle(paste0(varname_inplot, ", q1(t), factual = historical + rcp85, counterfactual = historicalNat")) +
  ylab("q1")
plot(p)
dev.off()

CvM_paris_filtered <- keep_onemodel_perinstitute(p12_allchain$CvM_df, unique(pr_cmip5[, c("institute", "model")]))
CvM_paris_filtered <- subset(CvM_paris_filtered, model != "HadCRUT")
CvM_paris_filtered$weight   <- CvM_paris_filtered$weight /  sum(CvM_paris_filtered$weight)
CvM_paris_filtered <- within(
  CvM_paris_filtered,
  model <- factor(model,  levels=model[order(weight, decreasing = TRUE)])
)
CvM_paris_filtered <- melt(CvM_paris_filtered, id.vars = c("institute", "model"))
pdf(paste0(cmip5_prefix,"_oxford_CvM.pdf"), width = 20/2.54, height = 20/2.54)
p <- ggplot(data = CvM_paris_filtered, aes(x = model, y = value, fill = grepl("(CM6)", model))) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(variable),  scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("darkgrey", "red"), labels = c("CMIP5", "CMIP6"), name= "CMIP")
plot(p)
dev.off()

pr_ghcnd_counterfactual_paris <- subset(pr_cmip5, model == "HadCRUT" & year <= 1900 & year > 1850, pr) %>% unlist()
pdf(paste0(cmip5_prefix,"_oxford_hist_checkH0.pdf"), width = 20/2.54, height = 20/2.54)
p <- hist_checkH0(
  X = pr_ghcnd_counterfactual_paris,
  lUhat = p12_allchain$lUhat[lmodels != "HadCRUT"],
  lmodels = lmodels[lmodels != "HadCRUT"],
  linstitutes = linstitutes[linstitutes != "HadCRUT"]
) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
plot(p)
dev.off()

pdf(paste0(cmip5_prefix,"_oxford_qqplot_checkH0.pdf"), width = 20/2.54, height = 25/2.54)
p <- qqplot_checkH0(
  X = pr_ghcnd_counterfactual_paris,
  lUhat = p12_allchain$lUhat[lmodels != "HadCRUT"],
  lmodels = lmodels[lmodels != "HadCRUT"],
  linstitutes = linstitutes[linstitutes != "HadCRUT"]
) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
plot(p)
dev.off()
# # Question de philippes sur les dates des max de Uhat
# # mapply(function(indice, model){subset(pr_cmip5, model == model & year <= 1900)[indice, "year"]}, indice = sapply(lp12_paris$lUhat[lmodels != "ghcnd"], which.max), model = lmodels[lmodels != "ghcnd"])
# gc(TRUE)

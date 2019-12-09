library(Rsolnp)

# solnp(c(5,5,5),
#       +       opt_func,
#       +       eqfun=equal,
#       +       eqB=15,
#       +       LB=c(0,0,0),
#       +       UB=c(100,100,100))
KLgauss <- function(mu1, mu2, sig21, sig22){
  -1/2 * (log(sig21) - log(sig22)) + (sig21 + (mu1 - mu2)^2) / (2 * sig22) - 1/2
}
gen_KLfun <- function(q1, V, I){
  stopifnot(length(q1) == length(V) & length(q1) == length(I))
  M <- length(q1)
  sig2_ref <- 0.5 * sum(1/I * 1/M^2)
  mu_ref <- 0.5
  KLfunc <- function(weight){
    # -1/2 * (log(sum(weight^2 * V / I)) - log(sig2_ref)) + (sum(weight^2 * V / I) + (sum(weight * q1) - mu_ref)^2) / (2 * sig2_ref) - 1/2  
    # + sum(weight^2 * I * KLgauss(q1, 0.5, V, 0.5))
    (sum(weight * q1) - mu_ref)^2 + sum(weight^2 * V / I) + sig2_ref - 2 * sqrt(sum(weight^2 * V / I) * sig2_ref)
  }
}

KLfunc <- gen_KLfun(q1 = rep(0.5, 2), V = rep(0.5, 2), I = rep(1, 2))
KLfunc(c(0.5, 0.5))

iter_KLfixpoint <- function(weights, q1, V, I){
  M <- length(q1)
  sig2_ref <- 0.5 * sum(1/I * 1/M^2)
  mu_ref <- 0.5
  mu <- sum(weights * q1)
  sig2 <- sum(weights^2 * V / I)
  sum_iV <- sum(1 / (V/I))
  (mu - mu_ref) / (V/I) * (sig2 / (sig2_ref - sig2)) * (q1 - (1/sum_iV) * sum(q1 / (V/I))) + 1/((V/I) * sum_iV)
}


M = 2
solnp(
  pars = c(0.2, 0.8), 
  fun = KLfunc,
  eqfun = function(weight) sum(weight),
  eqB = 1,
  LB=rep(0, M),
  UB=rep(1, M)
)

estimate_q1_stationnary <- function(z, x){
  q1 <- mean(ecdf(x)(z))
  return(q1)
}
estimate_q1max_stationnary <- function(z, x){
  q1max <- 0
  for(i in seq_along(x)){
    for(j in seq_along(x)){
      q1max <- q1max +  sum(max(x[i], x[j]) < z)
    }
  }
  q1max <- q1max / (length(x)^2 * length(z))
  return(q1max)
}

estimate_Vm_stationnary <- function(z, x){
 q1 <- estimate_q1_stationnary(z = z, x = x)
 q1max <- estimate_q1max_stationnary(z = z, x = x)
 expectation <- 1 - estimate_q1max_stationnary(z = x, x = z)
 Vm = q1max + expectation - 2 * q1^2
 return(Vm)
}

prepare_foroptim <- function(z, x){
  c("q1m" = estimate_q1_stationnary(z = z, x = x),
    "Vm" = estimate_Vm_stationnary(z = z, x = x),
    "Im" = length(z) * length(x)
  )
}




compute_allchain_KLweight <- function(lmodels, tas_cmip5, kernel, bandwidth, yearmin = 1850, yearmax = 1900){
  lp12 <- lapply(lmodels, function(model, tas_cmip5, kernel, bandwidth){
    print("---------------------------")
    print(model)
    tas_model <- tas_cmip5[tas_cmip5$model == model, ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical" | tas_model$experiment == "rcp85", ]
    print(head(tas_model_factual))
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", ]
    print(head(tas_model_counterfactual))
    # tas_model_counterfactual <- tas_model[tas_model$experiment == "historical" & tas_model$year <= 1900, ]
    tpred <- min(tas_model$year):max(tas_model$year)
    p12 <- estim_p12_ns(x = tas_model_counterfactual$tas, t = tas_model_factual$year, z = tas_model_factual$tas, tpred = tpred, kernel = kernel, h = bandwidth)
  }, tas_cmip5 = tas_cmip5, kernel = kernel, bandwidth = bandwidth)
  
  optimpars <- sapply(lmodels, function(model, tas_cmip5, yearmin, yearmax){
    # browser()
    print("---------------------------")
    # print(model)
    tas_model <- tas_cmip5[tas_cmip5$model == model & tas_cmip5$year >= yearmin & tas_cmip5$year <= yearmax , ]
    tas_model_factual <- tas_model[tas_model$experiment == "historical", "tas"]
    tas_model_counterfactual <- tas_model[tas_model$experiment == "historicalNat", "tas"]
    foroptim <- prepare_foroptim(z = tas_model_factual, x = tas_model_counterfactual)
    return(foroptim)
  }, tas_cmip5 = tas_cmip5, yearmin = 1850, yearmax = 1900)
  print(optimpars)
  # browser()
  KLfunc <- gen_KLfun(q1 = optimpars["q1m", ], V = optimpars["Vm", ], I = optimpars["Im", ])
  M <- length(lmodels)
  pars <- 1 / KLgauss(optimpars["q1m", ], 0.5, optimpars["Vm", ], 0.5)
  pars <- pars / sum(pars)
  pars <- rep(1/M, M)
  KLoptim <- solnp(
    pars = pars, 
    fun = KLfunc,
    eqfun = function(weight) sum(weight),
    eqB = 1,
    LB = rep(0, M),
    UB = rep(1, M)
  )

  p12_multimodel <- multimodel_average(lp12, lmodels, lweight = KLoptim$pars)  
  
  list(
    lp12 = lp12,
    KLoptim = KLoptim,
    lweight = KLoptim$pars,
    optimpars = optimpars,
    p12_multimodel = p12_multimodel
  )
}

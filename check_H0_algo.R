library(magrittr)

ecdf_nw <- function(q, tpred, t, x, kernel, bandwidth){
  dvect <- abs(tpred - t)
  kvect <- kernel(dvect, h = bandwidth)
  indicator <- (x < q)
  weighted.mean(indicator, w = kvect)
}

quantile_nw <- function(p, tpred, t, x, kernel, bandwidth){
  # cat("p = ", p , "; tpred = ", tpred,  " \n")
  dvect <- abs(tpred - t)
  kvect <- kernel(dvect, h = bandwidth)
  uniroot(f = function(q) {
    indicator <- (x < q)
    Gq <- weighted.mean(indicator, w = kvect)
    Gq - p
  }, interval = range(x) + c(-sd(x), +sd(x)))$root
}

compute_Uhat <- function(t, tm, Xm, Z, Zm, kernel, bandwidth){
  FmXm <- sapply(seq_along(tm), function(i){
    ecdf_nw(Xm[i], tm[i], tm, Zm, kernel = kernel_epanechnikov, bandwidth = bandwidth)
  })
  Uhat <- sapply(seq_along(tm), function(i){
    quantile_nw(FmXm[i], tm[i], t, Z, kernel = kernel_epanechnikov, bandwidth = bandwidth)
  })
  return(Uhat)
}

check_HO_qqplot <- function(X, Uhat){
  xyrange <- range(X, Uhat)
  qqplot(X, Uhat, main = "", xlim = xyrange, ylim = xyrange)
  abline(a = 0, b = 1, col = "red")
}

check_HO_rankhist <- function(X, Uhat){
  GmU <- ecdf(X)(Uhat)
  cut_GmU <- cut(GmU, seq(0, 1, 0.2), include.lowest = TRUE)
  # ks.test(GmU[GmU < 1 & GmU > 0], y = "punif")
  print(chisq.test(table(cut_GmU)))
  hist(GmU, breaks = 5, main = "")
  abline(h = length(GmU)/5, col = "red")
}

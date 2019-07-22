library(magrittr)
source("check_H0_algo.R")

n <- 2000
X <- rnorm(n)
Xm <- rnorm(n) * 2
# Z <- rnorm(n) + seq(0, 2, length.out = n)
# Zm <- rnorm(n) + seq(0, 2, length.out = n) + 2 
Z <- rnorm(n) + 1
Zm <- (rnorm(n) + 1)  * 3

x <- seq(min(X, Xm, Z, Zm), max(X, Xm, Z, Zm), length.out = 1000)
plot(x, qnorm(pnorm(x, sd = 2), sd = 1), xlim = range(x), ylim = range(x))
points(x, qnorm(pnorm(x, sd = 2, mean = 2), mean = 1, sd = 1), col = "blue")
abline(a = 0, b = 1, col = "red")
abline(h = 0, col  = "red")
abline(v = 0, col  = "red")

matplot(cbind(X, Xm, Z, Zm), type = "l", lty = 1)

Uhat <- ecdf(Zm)(Xm) %>% quantile(Z, probs = .)

qqplot(X, Uhat)
abline(a = 0, b = 1)

hist(ecdf(X)(Uhat), breaks = 5)
GmU <- ecdf(X)(Uhat)
cut_GmU <- cut(GmU, seq(0, 1, 0.2), include.lowest = TRUE)
ks.test(GmU[GmU < 1 & GmU > 0], y = "punif")
chisq.test(table(cut_GmU))


# With kernel smoothing ---------------------------------------------------

Z <- rnorm(n) + seq(0, 1, length.out = n)
ecdf_nw(0, tpred = 1, t = seq_along(Z), x = Z, kernel = kernel_epanechnikov, bandwidth = 200)  

quantile_nw(0.75, tpred = 1, t = seq_along(Z), x = Z, kernel = kernel_epanechnikov, bandwidth = 30)  

n <- 2000
t <- seq.int(n)
X <- rnorm(n)
Xm <- rnorm(n) * 2
Z <- rnorm(n) + seq(0, 2, length.out = n)
Zm <- (rnorm(n) + seq(0, 2, length.out = n)) * 2
# Z <- rnorm(n) + 1
# Zm <- (rnorm(n) + 1) * 2  

Uhat <- compute_Uhat(t, Xm, Z, Zm, kernel = kernel_epanechnikov, bandwidth = 200)
qqplot(X, Uhat)
abline(a = 0, b = 1)

check_HO_qqplot(X, Uhat)

check_HO_rankhist(X, Uhat)


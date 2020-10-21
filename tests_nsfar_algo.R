### Gumble case
library(magrittr)

kernel_gauss <- function(x, h = 1){
  dnorm(x / h) / h
}
attr(kernel_gauss, "K2_integrated") <-  1 / (2 * sqrt(pi))

kernel_epanechnikov <- function(x, h = 1){
  pmax(3/4 * (1 - (x/h)^2) / h, 0)
}
attr(kernel_epanechnikov, "K2_integrated") <-  3/5

estim_p12_ns <- function(x, t, z, tpred = sort(unique(t)), kernel, h = length(t)^(-1/5)){
  n <- length(t)
  m <- length(x)
  Gm <- ecdf(x)
  GmZ <- Gm(z)
  p12 <- sapply(tpred, function(ti){
    dvect <- abs(ti - t)
    kvect <- kernel(dvect, h = h)
    weighted.mean(GmZ, kvect)
  })
  p12t <- p12[match(t, table = tpred)]
  sig2 <- sapply(tpred, function(ti){
    dvect <- abs(ti - t)
    kvect <- kernel(dvect, h = h)
    weighted.mean((p12t - GmZ)^2, kvect)
  })
  if( !is.null(attr(kernel, "K2_integrated" ))) {
    K22 <- attr(kernel, "K2_integrated")
  } else{
    K22 <- integrate(function(x) kernel(x)^2, lower = -Inf, upper = Inf)$value
  }
  ft <- sapply(tpred, function(ti){
    dvect <- abs(ti - t)
    kvect <- kernel(dvect, h = h)
    mean(kvect)
  })
  p12_var <- 1/(n * h) * sig2/ft * K22
  # p12_var <-  p12_var + 1/(n * m * h) * p12 * (1 - p12)
  ## Part of Brownian Bridge
  minGmZ <-  outer(GmZ, GmZ, pmin)
  prodGmZ <- outer(GmZ, GmZ, "*")
  BGmZ_var <- sapply(tpred, function(ti){
    dvect <- abs(ti - t)
    kvect <- kernel(dvect, h = h)
    prodKj <- outer(kvect/sum(kvect), kvect/sum(kvect), "*")
    sum(prodKj * (minGmZ - prodGmZ))
  })
  p12_var <- p12_var + BGmZ_var / m
  list(
    tpred = tpred,
    p12_hat = p12,
    sigma_p12_hat = sqrt(p12_var),
    GmZ = GmZ
  )
}

estim_theta_ns <- function(p12_hat, sigma_p12_hat) {
  theta_hat <- 1/p12_hat - 1
  var_theta_hat <- sigma_p12_hat^2 / p12_hat^4
  list(
    theta_hat = theta_hat,
    sigma_theta_hat = sqrt(var_theta_hat)
  )
}

estim_theta.nswexp <- function(x, t, z, tpred = sort(unique(t)), kernel, h = NULL){
  #
  #  x (counterfactual), z (factual), rp = return periods

  # rounding numbers to be added to discrete values of x and z
  # rounding=seq(from = -0.499, to= 0.499, length = 2 *  max(length(z), length(x)))
  # z <- jitter_dat(z, rounding)
  # x <- jitter_dat(x, rounding)
  if (is.null(h)) {
    # Leave-One-Out Cross-Validation
    # h <- optimize(cv_loo_h,
    #               interval = c(
    #                 min(diff(sort(t))),
    #                 abs(diff(range(t)))
    #               ),
    #               x = x, t = t, z = z,
    #               kernel = kernel)$minimum
    h_totest <- seq(
      min(diff(tpred)),
      abs(diff(range(t))) / 4,
      by = median(diff(tpred))
    )
    cv_score <- sapply(
      h_totest,
      cv_loo_h, x = x, t = t, z = z, kernel = kernel
    )
    h <- h_totest[which.min(cv_score)]
    # h <- choose_h_for_wexp(x = x, t = t, z = z, kernel = kernel)$minimum
  }
  p12 <- estim_p12_ns(x = x, t = t, z = z, tpred = tpred, kernel = kernel, h = h)
  theta <- estim_theta_ns(p12_hat = p12$p12_hat, sigma_p12_hat = p12$sigma_p12_hat)

  # Creating W = - log(Ghat(Z)) with Ghat(Z)= average(K(Z-X_i)) for the kernel based on the arctan
  # GZhat <- estim_GZ.kernel(x, z)
  # GZhat <- estim_GZ.empirical(x, z)
  # W <- estim_W(GZhat)
  # W_normalized <- normalize_W(W, theta$theta_hat)
  # W_normalized <- replace_zeros(W_normalized)
  # Cox and Oakes Test for the exponential distribution
  # co_test <- farr::CoxOakes(x = W_normalized)
  theta_with_duplicate <- merge(
    data.frame(t = t, z = z),
    data.frame(t = p12$tpred, theta_hat = theta$theta_hat),
    by = "t"
  )
  utest <- kstest_unif(x = x, z = theta_with_duplicate$z, theta = theta_with_duplicate$theta_hat)
  names(theta$theta_hat) <- paste("theta_t", p12$t_unique, sep = "_")
  theta_fit <- list(
    x = x, t = t, z = z,
    tpred = p12$tpred,
    theta_hat = theta$theta_hat,
    sigma_theta_hat = theta$sigma_theta_hat,
    GmZ = p12$GmZ,
    utest_pvalue = utest$p.value,
    kernel = kernel, h = h
  )
  class(theta_fit) <- c("thetafitns_wexp")
  return(theta_fit)
}

coef.thetafitns_wexp <- function(object, ...){
  coef <- object$theta_hat
  names(coef) <- names(object$theta_hat)
  return(coef)
}

vcov.thetafitns_wexp <- function(object, ...){
  if(length(object$sigma_theta_hat) == 1){
    vcov <- matrix(object$sigma_theta_hat, nrow = 1, ncol = 1)
  } else{
    vcov <- diag(object$sigma_theta_hat^2)
  }
  rownames(vcov) <- colnames(vcov) <- names(object$theta)
  return(vcov)
}

hist.thetafitns_wexp <- function(x, ...){
  #
  # Checking (histogram) if -log(G(Z)) follows an exponentialdistribution
  #
  GmZ <- x$GmZ
  theta_df <- merge(
    data.frame(t = x$t, GmZ = x$GmZ),
    data.frame(t = x$tpred, theta_hat = x$theta_hat),
    by = "t"
  )
  theta_df <- theta_df[!(theta_df$GmZ <= 0 | theta_df$GmZ >= 1), ]
  GmZ_normalized <- theta_df$GmZ^(1/theta_df$theta_hat)
  hist(GmZ_normalized,
       freq = F,
       xlab = "Gm(Z)^(1/theta)",
       xlim =  range(GmZ_normalized),
       main = "Uniform fit for GmZ_normalized",
       breaks = max(9, length(GmZ_normalized)/10),
       col = "lightgray"
  )
  # grid(lwd = 3)
  box()
  xx <- seq(from = 0, to = max(GmZ_normalized), length = 200)
  lines(xx, dunif(xx, min = 0, max = 1), col = "darkblue", lwd = 3, lty = 2)
}

qqplot.thetafitns_wexp <- function(x, ...){
  #
  # Checking (qqplot) if -log(G(Z)) follows an exponentialdistribution
  #
  GmZ <- x$GmZ
  ci95 <- confint(x)
  ci95 <- merge(
    data.frame(t = x$t, GmZ = x$GmZ),
    data.frame(t = x$tpred, theta_hat = x$theta_hat, q025 = ci95[, 1], q975 = ci95[, 2]),
    by = "t"
  )
  ci95 <- ci95[!(ci95$GmZ <= 0 | ci95$GmZ >= 1), ]
  GmZ_normalized <- ci95$GmZ^(1/ci95$theta_hat)
  GmZ_normalized_inf <- ci95$GmZ^(1/ci95$q025)
  GmZ_normalized_sup <- ci95$GmZ^(1/ci95$q975)
  ll <- length(GmZ_normalized)
  pp <- ((1:ll) - 0.5) / ll
  expected <- qunif(pp, min = 0, max = 1)
  observed = sort(GmZ_normalized)
  observed_sup = sort(GmZ_normalized_sup)
  observed_inf = sort(GmZ_normalized_inf)
  xylim <- range(
    expected[is.finite(expected)],
    observed[is.finite(observed)]
  )
  plot(observed, expected,
       xlim = xylim, ylim = xylim,
       xlab = "Observed", ylab = "Expected",
       main = "Uniform QQ plot for \n GmZ_normalized=Gm(Z)^(1/theta)",
       col = "darkblue", pch = 20, cex = 1)
  # grid(lwd = 3)
  box()
  polygon(
    c(observed_inf, observed_sup[ll:1]),
    c(expected, expected[ll:1]),
    col = rgb(red = 0,green = .0, blue = .5, alpha = .3),
    border = F
  )
  abline(a = 0, b = 1 , col = "red", lwd = 3)
  points(observed, expected, col = "darkblue", pch = 20, cex = 1)
}

ecdf.thetafitns_wexp <- function(x, ...){
  GmZ <- x$GmZ
  theta_df <- merge(
    data.frame(t = x$t, GmZ = x$GmZ),
    data.frame(t = x$tpred, theta_hat = x$theta_hat),
    by = "t"
  )
  theta_df <- theta_df[!(theta_df$GmZ <= 0 | theta_df$GmZ >= 1), ]
  GmZ_normalized <- theta_df$GmZ^(1/theta_df$theta_hat)
  xlim <- range(GmZ_normalized)
  plot(ecdf(GmZ_normalized),
       xlim = xlim,
       xlab = "GmZ_normalized=Gm(Z)^(1/theta)",
       ylab = "",
       main = "Theoretical and Empirical CDFs", ...)
  # grid(lwd = 3);box()
  xx <- seq(from = 0 , to = xlim[2], length = 200)
  lines(xx, punif(xx, min = 0, max = 1), col = "darkgray", lwd = 4, lty = 2)
}

plot.thetafitns_wexp <- function(x, ...){
  order_t <- order(x$tpred)
  t_ordered <- x$t[order_t]
  theta_ordered <- x$theta_hat[order_t]
  theta_ci_ordered <- confint(x)[order_t,]
  plot(x$t, 1/x$GmZ - 1, pch = 20, col = "lightgrey",
       ylab = "theta_hat", xlab = "t", main = "theta(t)")
  lines(t_ordered, theta_ordered, lty = 2)
  matlines(t_ordered, theta_ci_ordered, col = "black",  lty = 1)
}

estim_farr_ns <- function(theta_hat, sigma_theta_hat, rp) {
  stopifnot(length(rp) == 1)
  farr_hat <- (1 - theta_hat) * (1 - 1/rp)
  farr_hat_var <- (1 - 1/rp)^2 * sigma_theta_hat^2
  cbind(farr_hat = farr_hat, sigma_farr_hat = sqrt(farr_hat_var))
}

estim_farr.nswexp <- function(object, rp){
  t <- object$t
  theta_hat <- object$theta_hat
  sigma_theta_hat <- object$sigma_theta_hat
  farr_arr <- vapply(rp,
                    FUN = estim_farr_ns,
                    FUN.VALUE = matrix(1, nrow = length(theta_hat), ncol = 2),
                    theta_hat = theta_hat, sigma_theta_hat = sigma_theta_hat
  )
  dimnames(farr_arr) <- list(t = t,
                            estim = c("far_hat", "sigma_far_hat"),
                            rp = rp)
  farr_hat <- matrix(farr_arr[, 1,], nrow = length(t), ncol = length(rp))
  sigma_farr_hat <- matrix(farr_arr[, 2,], nrow = length(t), ncol = length(rp))
  farr_fit <- list(t = t,
                   rp = rp,
                   farr_hat = farr_hat,
                   sigma_farr_hat = sigma_farr_hat,
                   GmZ = object$GmZ)
  class(farr_fit) <- c("farrfitns_wexp")
  return(farr_fit)
}

extract_rp <- function(object, rp){
  UseMethod("extract_rp", object)
}

extract_rp.farrfitns_wexp <- function(object, rp){
  irp <- which(object$rp %in% rp)
  farr_fit <- list(t = object$t,
                   rp = object$rp[irp],
                   farr_hat = object$farr_hat[, irp, drop = FALSE],
                   sigma_farr_hat = object$sigma_farr_hat[, irp, drop = FALSE],
                   GmZ = object$GmZ)
  class(farr_fit) <- c("farrfitns_wexp")
  return(farr_fit)
}

extract_rp.boot_farrfitns_wexp <- function(object, rp){
  irp <- which(object$rp %in% rp)
  farr_fit <- list(t = object$t,
                   rp = object$rp[irp],
                   farr_boot = object$farr_boot[, irp,, drop = FALSE],
                   theta_boot = object$theta_boot,
                   xboot = object$xboot,
                   tboot = object$tboot,
                   zboot = object$zboot)
  class(farr_fit) <- class(object)
  return(farr_fit)
}

coef.farrfitns_wexp <- function(object, ...){
  coef <- object$farr_hat
  dimnames(coef) <- list(t = object$t, rp = object$rp)
  return(coef)
}

vcov.farrfitns_wexp <- function(object, ...){
  vcov <- object$sigma_farr_hat^2
  dimnames(vcov) <- list(t = object$t, rp = object$rp)
  return(vcov)
}

confint.farrfitns_wexp <- function(object, parm, level = 0.95,  ...){
  alpha <- (1 - level) /2
  mu <- coef(object)
  sigma <- sqrt(vcov(object))
  if(!missing(parm)){
    irp <- seq_along(object$rp)
    it <- seq_along(object$t)
    if(!is.null(parm$rp)){
      irp <- which(object$rp %in% parm$rp)
    }
    if(!is.null(parm$it)){
      it <- which(object$rp %in% parm$t)
    }
    mu <- mu[it, irp]
    sigma <- sigma[it, irp]
  }
  ci_sup <- qnorm(mean = mu, sd = sigma, p = level + alpha)
  ci_inf <- qnorm(mean = mu, sd = sigma, p = alpha)
  ci_arr <- array(c(ci_inf, ci_sup), dim = c(dim(ci_sup), 2))
  dimnames(ci_arr) <- list(t = rownames(mu),
                           rp = colnames(mu),
                           ci = c(alpha, level + alpha))
  return(ci_arr)
}
get_rp.farrfitns_wexp <- function(object, ...){
  object$rp
}

get_t <- function(object, ...){
  UseMethod("get_t", object)
}

get_t.farrfitns_wexp <- function(object, ...){
  object$t
}

get_t.boot_farrfitns_wexp <- function(object, ...){
  object$tboot[[1]]
}

get_farr.farrfitns_wexp <- function(object, ...){
  object$farr_hat
}

get_farr.boot_farrfitns_wexp <- function(object, ...){
  apply(object$farr_boot, 1:2, mean)
}

get_GmZ <- function(object, ...){
  UseMethod("get_GmZ", object)
}

get_GmZ.farrfitns_wexp <- function(object, ...){
  object$GmZ
}

get_GmZ.boot_farrfitns_wexp <- function(object, ...){
  ecdf(object$xboot[[1]])(object$zboot[[1]])
}

plot.farrfitns_wexp <- function(x, plot_ci = TRUE, ...){
  rp <- get_rp(x)
  order_t <- get_t(x) %>% order()
  t_ordered <- get_t(x)[order_t]
  t_unique <- unique(t_ordered)
  it_unique  <- vapply(t_unique,
                       FUN = function(x) min(which(x == t_ordered)),
                       FUN.VALUE = 1)
  farr <- get_farr(x)[order_t,, drop = FALSE]
  farr <- farr[it_unique,, drop = FALSE]
  farr_ci_ordered <- confint(x)[order_t,,, drop = FALSE]
  farr_ci_ordered <- farr_ci_ordered[it_unique,,, drop = FALSE]
  tmat <- matrix(rep(t_ordered, length(rp)),
                 nrow = length(t_ordered),
                 ncol = length(rp))
  rpmat <- matrix(rep(rp, length(t_ordered)),
                  nrow = length(t_ordered),
                  ncol = length(rp),
                  byrow = TRUE)
  tmat_unique <- matrix(rep(t_unique, length(rp)),
                 nrow = length(t_unique),
                 ncol = length(rp))
  rpmat_unique <- matrix(rep(rp, length(t_unique)),
                  nrow = length(t_unique),
                  ncol = length(rp),
                  byrow = TRUE)
  rpmat_col <- cut(rpmat_unique, breaks = seq(min(rp)-0.1, max(rp + .1), length.out = length(rp) + 1))
  rp_col <- cut(rp, breaks = seq(min(rp)-0.1, max(rp + .1), length.out = length(rp) + 1))
  itime <- sapply(seq.int(length(t_unique)-1),
                  function(i){
                    c(i, i + 1, i +1, i)
                  })
  t_quads3d <- t_unique[itime]
  ciband_byrp <- function(irp){
    farr_quads3d <- sapply(seq.int(length(t_unique)-1),
                           function(i){
                             c(farr_ci_ordered[i, irp, 1],
                               farr_ci_ordered[i + 1, irp, 1],
                               farr_ci_ordered[i + 1, irp, 2],
                               farr_ci_ordered[i, irp, 2])
                           })
    rgl::quads3d(x = t_quads3d,
                 y = rep(rp[irp], 2 * length(t_unique)),
                 z = farr_quads3d,
                 col = rainbow(length(rp))[rp_col[irp]], alpha = 0.3)
  }
  rgl::plot3d(tmat_unique, rpmat_unique, farr,
              col = rainbow(length(rp))[rpmat_col],
              xlab = "t", ylab = "rp", zlab = "far")
  if(plot_ci){
    sapply(seq_along(rp), ciband_byrp)
  }
  if(length(rp) == 1){
    farr_points <- (1 - (1/x$GmZ - 1)) * (1 - 1/max(rpmat))
    farr_points[is.infinite(farr_points)] <- NA
    #     farr_points[farr_points <  min(farr_ci_ordered)  * 1.2] <- NA
    #     farr_points[farr_points >  max(farr_ci_ordered) * 1.2] <- NA
    farr_points[farr_points <  -0.5] <- NA
    farr_points[farr_points >  1.2] <- NA
    rgl::points3d(tmat, rpmat, farr_points)
    rgl::grid3d(c("y+"))
    rgl::aspect3d(1, 1, 0.75)
  } else{
    rgl::grid3d(c("x", "y+", "z"))
  }
}

plot.farrfitns_wexp <- function(x, plot_ci = TRUE, ...){
  ellipsis <- list(...)
  if(is.null(ellipsis$level)){
    ellipsis$level  <-  0.95
  }
  rp <- get_rp(x)
  order_t <- get_t(x) %>% order()
  t_ordered <- get_t(x)[order_t]
  t_unique <- unique(t_ordered)
  it_unique  <- vapply(t_unique,
                       FUN = function(x) min(which(x == t_ordered)),
                       FUN.VALUE = 1)
  farr <- get_farr(x)[order_t,, drop = FALSE]
  farr <- farr[it_unique,, drop = FALSE]
  farr_ci_ordered <- confint(x, level = ellipsis$level)[order_t,,, drop = FALSE]
  farr_ci_ordered <- farr_ci_ordered[it_unique,,, drop = FALSE]
  if(is.null(ellipsis$ylim)){
    ellipsis$ylim <- range(farr,  farr_ci_ordered, na.rm = TRUE)
  }
  tmat <- matrix(rep(t_ordered, length(rp)),
                 nrow = length(t_ordered),
                 ncol = length(rp))
  rpmat <- matrix(rep(rp, length(t_ordered)),
                  nrow = length(t_ordered),
                  ncol = length(rp),
                  byrow = TRUE)
  tmat_unique <- matrix(rep(t_unique, length(rp)),
                 nrow = length(t_unique),
                 ncol = length(rp))
  rpmat_unique <- matrix(rep(rp, length(t_unique)),
                  nrow = length(t_unique),
                  ncol = length(rp),
                  byrow = TRUE)
  rp_col <- cut(rp, breaks = seq(min(rp) - 0.1, max(rp + .1), length.out = length(rp) + 1))
  ciband_byrp <- function(irp){
    x <- c(t_unique, rev(t_unique))
    y <- c(farr_ci_ordered[, irp, 1],
           rev(farr_ci_ordered[, irp, 2]))
    col <- rainbow(length(rp))[rp_col[irp]] %>%
      # makeTransparent(alpha = 50)
      adjustcolor(alpha.f=0.2)
    polygon(x = x,  y = y, col = col, border = FALSE)
  }
  plot_onerp <- function(irp){
    if (is.null(ellipsis$main)) {
      ellipsis$main  <-  "far_rp(t)"
    }
    plot(tmat_unique[, irp], farr[, irp],
          type = "l", lwd = 2, lty = 1,
          ylim = ellipsis$ylim,
          col = rainbow(length(rp))[rp_col[irp]],
          xlab = "t", ylab = "far",
          main = ellipsis$main,
          sub = paste("rp =", rp[irp], ", CI level =", ellipsis$level))
    ciband_byrp(irp)
    farr_points <- (1 - (1/get_GmZ(x) - 1)) * (1 - 1/rp[irp])
    points(tmat[, irp], farr_points)
    grid()
  }
  compute_nrowncol <- function(rp){
    lrp <- length(rp)
    nrow <- floor(sqrt(length(rp)))
    ncol <- ceiling(lrp / nrow)
    c(nrow, ncol)
  }
  if(length(rp) > 1){
    par(mfrow = compute_nrowncol(rp))
  }
  sapply(seq_along(rp),
         function(irp){
           plot_onerp(irp)
         })
}

#note: always pass alpha on the 0-255 scale
makeTransparent <- function(someColor, alpha=100){
  newColor <- col2rgb(someColor)
  apply(newColor, 2,
        function(curcoldata){
          rgb(red = curcoldata[1],
              green = curcoldata[2],
              blue = curcoldata[3],
              alpha = alpha,
              maxColorValue=255)
        })
}

# Bootstrap estimation
boot_p12 <- function(x, t, z, kernel = kernel_epanechnikov, h = length(t)^(-1/5), B = 100){
  xboot <- lapply(seq.int(B), function(i) sample(x = x, replace = TRUE))
  xboot[[1]] <- x
  n <- length(t)
  t_unique <- sort(unique(t))
  n_unique <- length(t_unique)
  itboot <- lapply(seq.int(B), function(i){
    sample.int(n,
               size = n,
               replace = TRUE)
  })
  tboot <- lapply(itboot, function(it) t[it])
  tboot[[1]] <- t
  zboot <- lapply(itboot, function(it) z[it])
  zboot[[1]] <- z
  p12_boot <- mapply(function(x, t, z){
    dmat <- c(t, t_unique) %>% matrix(ncol = 1) %>%
      dist(method = "euclidian") %>%
      as.matrix()
    dmat <- dmat[1:n, (n+1):(n + n_unique)]
    kmat <- apply(dmat, 2, kernel, h = h)
    Gm <- ecdf(x)
    GmZ <- Gm(z)
    p12 <- apply(kmat, 2, function(weight){
      weighted.mean(GmZ, weight)
    })
    return(p12)
  }, x = xboot, t = tboot, z = zboot)
  list(
    t_unique = t_unique,
    p12_boot = p12_boot,
    xboot = xboot,
    tboot = tboot,
    zboot = zboot
  )
}

boot_theta <- function(p12_boot){
  1/p12_boot - 1
}

boot_theta_fit.nswexp <- function(x, t, z,
                                kernel = kernel_epanechnikov,
                                h = length(t)^(-1/5),
                                B = 100){
  p12_boot <- boot_p12(x = x, t = t, z = z,
                       kernel = kernel, h = h,
                       B = B)
  theta_boot <- boot_theta(p12_boot$p12_boot)
  list(
    t_unique= p12_boot$t_unique,
    theta_boot = theta_boot,
    xboot = p12_boot$xboot,
    tboot = p12_boot$tboot,
    zboot = p12_boot$zboot
  )
}

boot_farr <- function(theta_boot, rp) {
  if(length(rp) != 1) stop("length != 1")
  (1 - theta_boot) * (1 - 1/rp)
}

boot_p1r <- function(theta_boot, rp){
  if(length(rp) != 1) stop("length != 1")
  1 / (1 + (rp - 1) * theta_boot)
}

create_boot_stat_fit.nswexp <- function(stat_name, boot_stat){
  boot_farr_fit.nswexp <- function(object , rp){
    stat_boot <- vapply(rp,
                        FUN = function(rpi){
                          boot_stat(rp = rpi, theta_boot = object$theta_boot)
                        },
                        FUN.VALUE = object$theta_boot)
    stat_boot <- aperm(stat_boot, c(1, 3, 2))
    dimnames(stat_boot) <- list(t = object$t_unique,
                                rp = rp,
                                B = seq.int(length(object$tboot)))
    stat_res <- list(stat_boot = stat_boot, rp = rp)
    names(stat_res) <- c(
      paste(stat_name, "_boot", sep = ""),
      "rp"
    )
    stat_res <- c(stat_res, object)
    class(stat_res) <- c(
      paste("boot_", stat_name, "fitns_wexp", sep = ""),
      paste(stat_name, "fitns_wexp", sep = "")
    )
    return(stat_res)
  }
}
boot_farr_fit.nswexp <- create_boot_stat_fit.nswexp(stat_name = "farr", boot_stat = boot_farr)
boot_p1r_fit.nswexp <- create_boot_stat_fit.nswexp(stat_name = "p1r", boot_stat = boot_p1r)

# boot_farr_fit.nswexp <- function(object , rp){
#   farr_boot <- vapply(rp,
#                       FUN = function(rpi){
#                         boot_farr(rp = rpi, theta_boot = object$theta_boot)
#                       },
#                       FUN.VALUE = object$theta_boot)
#   farr_boot <- aperm(farr_boot, c(1, 3, 2))
#   dimnames(farr_boot) <- list(t = object$tboot[[1]],
#                               rp = rp,
#                               B = seq.int(length(object$tboot)))
#   farr_res <- c(list(farr_boot = farr_boot, rp = rp), object)
#   class(farr_res) <- c("boot_farrfitns_wexp", "farrfitns_wexp")
#   return(farr_res)
# }

confint.boot_farrfitns_wexp <- function(object, parm, level = 0.95, ...){
  alpha <- (1 - level) / 2
  farr_boot <- object$farr_boot
  if(!missing(parm)){
    irp <- seq_along(object$rp)
    it <- seq_along(object$t)
    if(!is.null(parm$rp)){
      irp <- which(object$rp %in% parm$rp)
    }
    if(!is.null(parm$it)){
      it <- which(object$rp %in% parm$t)
    }
    farr_boot <- farr_boot[it, irp, ]
  }
  ci_sup <- apply(farr_boot, 1:2, quantile, probs = level + alpha)
  ci_inf <- apply(farr_boot, 1:2, quantile, probs = alpha)
  ci_arr <- array(c(ci_inf, ci_sup), dim = c(dim(ci_sup), 2))
  dimnames(ci_arr) <- list(t = object$tboot[[1]],
                           rp = object$rp,
                           ci = c(alpha, level + alpha))
  return(ci_arr)
}


################################################################################
# Test de distribution exponentielle de W

normalize_W <- function(W, theta){
  exp(log(W) - log(theta))
}

replace_zeros <- function(W_norm){
  izero <- which(W_norm == 0)
  if (length(izero) == 0) return(W_norm)
  minW <- min(W_norm[-izero])
  W_norm[izero] <- runif(length(izero), min = 0, max = minW)
  return(W_norm)
}

replace_ones <- function(Gz){
  ione <- which(Gz >= 1)
  if (length(ione) == 0) return(Gz)
  maxW <- max(Gz[-ione])
  Gz[ione] <- runif(length(ione), min = maxW, max = 1 - 1E-16)
  return(Gz)
}


# Bandwidth choice

choose_h_for_wexp <- function(x = x, t = t, z = z,
                              kernel = kernel_epanechnikov){
  GZhat <- estim_GZ.kernel(x, z)
  W <- estim_W(GZhat)
  ftomin <- function(h){
    p12 <- estim_p12_ns(x = x, t = t, z = z, kernel = kernel, h = h)
    theta <- with(p12, estim_theta_ns(p12_hat = p12_hat, sigma_p12_hat = sigma_p12_hat))
    W_normalized <- normalize_W(W, theta$theta_hat)
    W_normalized <- replace_zeros(W_normalized)
    pvalue <- farr::CoxOakes(x = W_normalized)$p.value
    return(1 - pvalue)
  }
  roughly <-  length(t)^(-1/5)
  optimize(ftomin, lower = 0, upper = max(t),
           tol = .Machine$double.eps^(1/100))
}

cv_loo_h <- function(x = x, t = t, z = z,
                   kernel = kernel_epanechnikov, h = h){
  Gm <- ecdf(x)
  GmZ <- Gm(z)
  ft <- sapply(t, function(ti){
    dvect <- abs(ti - t)
    kvect <- kernel(dvect, h = h)
    mean(kvect)
  })
  compute_res2_i <- function(i){
    dvect <- abs(t[i] - t[-i])
    kvect <- kernel(dvect, h = h)
    p12i <- weighted.mean(GmZ[-i], w = kvect)
    e2i <- (GmZ[i] - p12i)^2
    # e2i <- abs(GmZ[i] - p12i)
  }
  e2 <- vapply(seq_along(t),
               FUN = compute_res2_i,
               FUN.VALUE = numeric(1))
  # cvh <- weighted.mean(e2, w = ft)
  cvh <- weighted.mean(e2, w = (GmZ > 0 & GmZ < 1))
  return(cvh)
}



# Simulation of Randomly Weighted Brownian Bridge
drawBBNW <- function(tref, size, h){
  t <- runif(size)
  # t <- seq(0, 1, length.out = size + 2)[-c(1, size + 2)]
  mu = 2 + exp(t)
  # x = rgev(size * 5, loc = 0, scale = 1, shape = 0)
  z = rgev(size, loc = mu, scale = 1, shape = 0)
  # GmZ <- ecdf(x)(z)
  GZ <- pgev(z, loc = 0, scale = 1, shape = 0)
  GZ_sorted <- c(0, sort(GZ), 1)
  GZ_diff <- diff(GZ_sorted)
  W <- cumsum(rnorm(length(GZ_diff), sd = sqrt(GZ_diff)))
  W1 <- tail(W, 1)
  GZ_sorted <- GZ_sorted[-1]
  BGZ_sorted <- W - GZ_sorted * W1
  BGZ <- BGZ_sorted[-length(BGZ_sorted)]
  BGZ[order(GZ)] <- BGZ
  Kh <- kernel_gauss(sqrt((t-tref)^2), h = h)
  KhBGZ <- Kh * BGZ
  # plot(KhBGZ[order(t)])
  # plot(KhBGZ[-length(KhBGZ)], KhBGZ[-1])
  # scan(n =  1)
  # print(mean(KhBGZ == 0))
  # mean(KhBGZ)
  # KhBGZ[1:10]
  KhBGZ
}

varempirical <- function(size){
  h = size^(-1/5)/10
  samples_BBNW <- sapply(1:1000, function(i) drawBBNW(tref = tref, size = size, h = h))
  # sum(apply(samples_BBNW, 1, function(x) mean(x^2)))/size^2
  pairs(t(samples_BBNW[1:5, ]))
  cormat <- as.matrix(cor(t(samples_BBNW[1:1000, ])))
  hist(cormat[upper.tri(cormat)], breaks = seq(-1, 1, by = 0.01))
  mean(apply(samples_BBNW, 2,  mean)^2)
}

gen_h0 <- function(x, z, theta){
  x_h0 <- rgev(length(x), loc = 0, scale = 1, shape = 0)
  z_h0 <- rgev(length(z), loc = -log(theta), scale = 1, shape = 0)
  # z_h0 <- rgev(length(z), loc = mu, scale = 1, shape = 0)
  # GZhat_h0 <- estim_GZ.kernel(x_h0, z_h0)
  GZhat_h0 <- estim_GZ.empirical(x_h0, z_h0)
  W_h0 <- estim_W(GZhat_h0)
  # xfit_h0 <- fgev(x_h0)
  # W_h0 <- with(xfit_h0, -log(pgev(z_h0, loc = estimate[1], scale = estimate[2], shape = estimate[3])))
  # W_h0 <- -log(pgev(z_h0, loc = 0, scale = 1, shape = 0))
  W_h0_normalized <- normalize_W(W_h0, theta)
  # W_h0_normalized <- W_h0_normalized[W_h0_normalized > 0]
  ks.test(W_h0_normalized, "pexp", rate = 1)$statistic
}

compute_exptest_pvalue <- function(x, z, theta, n_h0 = 100){
  GZhat <- estim_GZ.empirical(x, z)
  W <- estim_W(GZhat)
  W_normalized <- normalize_W(W, theta)
  dist_ks <- ks.test(W_normalized, "pexp", rate = 1)$statistic
  dist_ks_h0 <- sapply(seq.int(n_h0), function(i) gen_h0(x, z, theta))
  # hist(dist_ks_h0, breaks = 50)
  # abline(v = dist_ks)
  pvalue <- mean(dist_ks < dist_ks_h0)
  return(pvalue)
}

## Tests that G(Z)^1/theta ~ Beta(1, 1) = U[0, 1]
## Equivalent to testing that -log(G(Z))/theta ~ Exp(1)
ptrunc_beta <- function(q, shape1, shape2, a, b, ncp = 0, lower.tail = TRUE, log.p = FALSE){
  ptrunc(q, spec = "beta", a = 0, b = 1,
         shape1 = shape1, shape2 = shape2,
         ncp = ncp, lower.tail = lower.tail, log.p = log.p)
}
dtrunc_beta <- function(x, shape1, shape2, a, b, ncp = 0, log = FALSE){
  dtrunc(x, spec = "beta", a = 0, b = 1,
         shape1 = shape1, shape2 = shape2,
         ncp = ncp, log = log)
}
kstest_unif <- function(x, z, theta){
  GZhat <- estim_GZ.empirical(x, z)
  # GZhat <- G_theo(z)
  GZhat_trunc <- GZhat[!(GZhat <= 0 | GZhat >= 1)]
  theta_trunc <- theta[!(GZhat <= 0 | GZhat >= 1)]
  # hist(GZhat_trunc^(1/theta_theo_trunc), freq = FALSE, breaks = 50)
  # xx <- seq(0, 1, 0.01)
  # lines(xx, dtrunc_beta(xx, shape1 = 1, shape2 = 1))
  GZhat_normalized <- GZhat_trunc^(1/theta_trunc)
  # ks.test(
  #   GZhat_trunc^(1/theta_trunc),"ptrunc_beta",
  #   shape1 = 1, shape2 = 1,
  #   a = min(GZhat_normalized),
  #   b = max(GZhat_normalized)
  # )
  ks.test(
    GZhat_normalized, "punif",
    min = min(GZhat_normalized),
    max = max(GZhat_normalized)
  )
}

cross_validate_h <- function(x, t, z, tpred, kernel){
  #
  #  x (counterfactual), z (factual), rp = return periods
    h_totest <- seq(
      min(diff(tpred)),
      abs(diff(range(t))) / 4,
      by = median(diff(tpred))
    )
    cv_score <- sapply(
      h_totest,
      cv_loo_h, x = x, t = t, z = z, kernel = kernel
    )
    h <- h_totest[which.min(cv_score)]
    # h <- choose_h_for_wexp(x = x, t = t, z = z, kernel = kernel)$minimum
    return(h)
}

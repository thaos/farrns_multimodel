library(magrittr)

qqplot_checkH0 <- function(lXmZm, lmodels, linstitutes) {
  adtest <- data.frame(
    institute = linstitutes,
    model = lmodels,
    pvalue = sapply(lXmZm, function(XmZm) twosamples::ad_test(XmZm$Xm, XmZm$Zm)[2])
  )
  qq_df <- mapply(
    function(XmZm, model, institute) {
      qqXmZm <- qqplot(XmZm$Xm, XmZm$Zm, plot.it = FALSE) %>%
        as.data.frame()
      names(qqXmZm) <- c("Xm", "Zm")
      cbind(institute = institute, model = model, qqXmZm)
    },
    XmZm = lXmZm, model = lmodels, institute = linstitutes,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  xylim <- range(qq_df$Xm, qq_df$Zm)
  ggplot(data = qq_df) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x = Xm, y = Zm)) +
    geom_point(
      data = adtest,
      aes(x = xylim[1] + 0.1 * diff(xylim),
          y = xylim[1] + 0.9 * diff(xylim),
          fill = pvalue, size = 2 * pvalue
      ),
      shape = 23
    ) +
  guides(size = "none") +
    facet_wrap(~ institute + model, ncol = ceiling(sqrt(length(lmodels)))) +
    scale_fill_gradientn(colours = rev(magma(5)), limits = c(0, 1)) +
    ggtitle("Checking A: qq-plot(Xm, Zm)") +
    coord_fixed(xlim = xylim, ylim = xylim)
}

dist_checkH0 <- function(lp12, lmodels) {
  dist_df <- mapply(
    function(p12, model) {
      p12 <- p12$p12_hat[lp12$tpred <= 1900 & p12$tpred >= 1900]
      d <- -length(p12) * log(2) - 1/2 * sum(log(p12 * (1-p12)))
      data.frame(model = model, dist = d) },
    p12 = lp12, model = lmodels,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)
  return(dist_df)
}


keep_onemodel_perinstitute <- function(dist_df, linstitutes) {
  dist_df <- cbind(institute = linstitutes, dist_df)
  for (institute in unique(linstitutes)) {
    i_institute <- dist_df$institute == institute
    mindist <- min(dist_df$dist[i_institute])
    dist_df$dist[i_institute & dist_df$dist > mindist] <- +Inf
  }
  return(dist_df)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

kldiv <- function(p12) return(-length(p12) * log(2) - 1/2 * sum(log(p12 * (1-p12))))
compute_kl_weights <- function(experts) {
  M <- ncol(experts)
  kloptim <- quiet(Rsolnp::solnp(
    pars = rep(1, M) / M, 
    fun = function(weight){
      kldiv(apply(experts, 1, weighted.mean, w = weight))
    },
    eqfun = function(weight) sum(weight),
    eqB = 1,
    LB = rep(0, M),
    UB = rep(1, M)
  ))
  return(kloptim$pars)
}

compute_oracle_weights <- function(experts) {
  ntimes <- nrow(experts)
  oracle <- opera::oracle(
    Y = rep(1 / 2, ntimes),
    experts = experts,
    model = "convex",
    loss.type = "square"
  )
  return(
    as.numeric(oracle$coefficients)
  )
}






multimodel_average <- function(lp12, lmodels, lweights = 1 / length(lmodels)) {
  p12_df <- mapply(function(p12, model, weight) {
    data.frame(
      model = model,
      weight = weight,
      year = p12$tpred,
      p12 = p12$p12_hat,
      sig = p12$sigma_p12_hat
    )
  }, p12 = lp12, model = lmodels, weight = lweights, SIMPLIFY = FALSE) %>%
    do.call(rbind, .)
  p12_byyear <- split(p12_df, f = p12_df$year)
  lapply(p12_byyear, function(df) {
    weight <- df$weight / sum(df$weight)
    p12 <- sum(df$p12 * weight)
    var_intra <- sum(df$sig^2 * weight^2)
    var_inter <- sum((df$p12 - p12)^2 * weight^2)
    sig <- sqrt(var_intra + var_inter)
    data.frame(year = unique(df$year), p12 = p12, sig = sig)
  }) %>% do.call(rbind, .)
}

indivEN.pois.exact <- function(z, n.eff, mean0, sd0, level=0.95,
                               grid.pi.null=seq(0.8,0.999,0.005),
                               grid.factor.mean.gamma=seq(-5,5,0.1),
                               grid.factor.var.gamma=seq(1,5,0.1), target=1) {
  fit.rlm <- MASS::rlm(z~1, method="M", scale.est="MAD", psi=MASS::psi.huber)
  if (missing(mean0)) mean0 <- unname(fit.rlm$coef)
  if (missing(sd0)) sd0 <- fit.rlm$s
  if (!is.numeric(mean0)) stop("Invalid mean0!")
  if (!is.numeric(sd0)) stop("Invalid sd0!")

  factor0.mean.gamma <- mean0/sqrt(median(n.eff))/sqrt(target)
  factor0.var.gamma <-
    uniroot(function(x) x + target*(x^4-x^2)*median(n.eff) - sd0^2,
            lower=0, upper=10)$root
  c.val <- qnorm((1-level)/2, lower.tail=F) # critical value
  mean.supp <- factor0.mean.gamma*sqrt(n.eff)*sqrt(target)
  sd.supp <- sqrt(factor0.var.gamma +
                    target*(factor0.var.gamma^4-factor0.var.gamma^2)*n.eff)
  lower0 <- mean.supp - c.val*sd.supp # lower limits
  upper0 <- mean.supp + c.val*sd.supp # upper limits
  idx.null <- which(z>=lower0 & z<=upper0)
  z.null <- z[idx.null]
  N.null <- length(z.null)

  negloglkd <- function(pi.null, factor.mean.gamma, factor.var.gamma) {
    means <- factor.mean.gamma * sqrt(n.eff) * sqrt(target)
    sds <- sqrt(factor.var.gamma +
                  target*(factor.var.gamma^4-factor.var.gamma^2)*n.eff)
    Q.nonnull <- pnorm((upper0[-idx.null]-means[-idx.null])/sds[-idx.null]) -
      pnorm((lower0[-idx.null]-means[-idx.null])/sds[-idx.null])
    negloglikelihood <-
      sum(log(sds[idx.null]) + (z.null - means[idx.null])^2/(2*(sds[idx.null])^2)) -
      N.null * log(pi.null) - sum(log(1-pi.null*Q.nonnull))
    return(ifelse(is.nan(negloglikelihood), Inf, negloglikelihood))
  }

  mat <- sapply(grid.pi.null, function(pi.null) {
    mat <- sapply(grid.factor.mean.gamma, function(factor.mean.gamma) {
      mat <- sapply(grid.factor.var.gamma, function(factor.var.gamma) {
        c(pi.null, factor.mean.gamma, factor.var.gamma,
          negloglkd(pi.null, factor.mean.gamma, factor.var.gamma))
      })
      mat[,which.min(tail(mat, n=1))]
    })
    mat[,which.min(tail(mat, n=1))]
  })
  opt.init <- mat[,which.min(tail(mat, n=1))]

  inc.mean.gamma <- min(diff(grid.factor.mean.gamma)) # incre of factor.mean.gamma seq
  inc.var.gamma <- min(diff(grid.factor.var.gamma)) # incre of factor.var.gamma
  grid.factor.mean.gamma <-
    seq(opt.init[2]-inc.mean.gamma,
        opt.init[2]+inc.mean.gamma, inc.mean.gamma/10)
  grid.factor.var.gamma <-
    seq(max(opt.init[3]-inc.var.gamma, 1),
        opt.init[3]+inc.var.gamma, inc.var.gamma/10)

  mat <- sapply(grid.factor.mean.gamma, function(factor.mean.gamma) {
    mat <- sapply(grid.factor.var.gamma, function(factor.var.gamma) {
      c(opt.init[1], factor.mean.gamma, factor.var.gamma,
        negloglkd(opt.init[1], factor.mean.gamma, factor.var.gamma))
    })
    mat[,which.min(tail(mat, n=1))]
  })
  opt.final <- mat[,which.min(tail(mat, n=1))]
  names(opt.final) <- c("pi.null", "factor.mean.gamma", "factor.var.gamma", "negloglkd.opt")
  return(opt.final)
}

indivEN.pois.approx <- function(z, n.eff, mean0, sd0, level=0.95,
                                grid.pi.null=seq(0.8,0.999,0.005),
                                grid.mean.z=seq(-5,5,0.1),
                                grid.var.gamma=seq(0,5,0.1), target=1) {

  fit.rlm <- MASS::rlm(z~1, method="M", scale.est="MAD", psi=MASS::psi.huber)
  if (missing(mean0)) mean0 <- unname(fit.rlm$coef)
  if (missing(sd0)) sd0 <- fit.rlm$s
  if (!is.numeric(mean0)) stop("Invalid mean0!")
  if (!is.numeric(sd0)) stop("Invalid sd0!")

  var0.gamma <- (sd0^2-1) / median(n.eff) / target
  c.val <- qnorm((1-level)/2, lower.tail=F) # critical value
  mean.supp <- mean0
  sd.supp <- sqrt(1 + var0.gamma*n.eff*target)
  lower0 <- mean.supp - c.val*sd.supp # lower limits
  upper0 <- mean.supp + c.val*sd.supp # upper limits
  idx.null <- which(z>=lower0 & z<=upper0)
  z.null <- z[idx.null]
  N.null <- length(z.null)

  negloglkd <- function(pi.null, mean.z, var.gamma) {
    sds <- sqrt(1 + var.gamma*n.eff*target)
    Q.nonnull <- pnorm((upper0[-idx.null]-mean.z)/sds[-idx.null]) -
      pnorm((lower0[-idx.null]-mean.z)/sds[-idx.null])
    negloglikelihood <-
      sum(log(sds[idx.null]) + (z.null - mean.z)^2/(2*(sds[idx.null])^2)) -
      N.null * log(pi.null) - sum(log(1-pi.null*Q.nonnull))
    return(ifelse(is.nan(negloglikelihood), Inf, negloglikelihood))
  }

  mat <- sapply(grid.pi.null, function(pi.null) {
    mat <- sapply(grid.mean.z, function(mean.z) {
      mat <- sapply(grid.var.gamma, function(var.gamma) {
        c(pi.null, mean.z, var.gamma,
          negloglkd(pi.null, mean.z, var.gamma))
      })
      mat[,which.min(tail(mat, n=1))]
    })
    mat[,which.min(tail(mat, n=1))]
  })
  opt.init <- mat[,which.min(tail(mat, n=1))]

  inc.mean.z <- min(diff(grid.mean.z)) # incre of mean.z seq
  inc.var.gamma <- min(diff(grid.var.gamma)) # incre of var.gamma
  grid.mean.z <-
    seq(opt.init[2]-inc.mean.z,
        opt.init[2]+inc.mean.z, inc.mean.z/10)
  grid.var.gamma <-
    seq(max(opt.init[3]-inc.var.gamma, 0),
        opt.init[3]+inc.var.gamma, inc.var.gamma/10)

  mat <- sapply(grid.mean.z, function(mean.z) {
    mat <- sapply(grid.var.gamma, function(var.gamma) {
      c(opt.init[1], mean.z, var.gamma,
        negloglkd(opt.init[1], mean.z, var.gamma))
    })
    mat[,which.min(tail(mat, n=1))]
  })
  opt.final <- mat[,which.min(tail(mat, n=1))]
  names(opt.final) <- c("pi.null", "mean.z", "var.gamma", "negloglkd.opt")
  return(opt.final)
}

indivEN.pois.0meanapprox <- function(z, n.eff, sd0, level=0.95,
                                     grid.pi.null=seq(0.8,0.999,0.005),
                                     grid.var.gamma=seq(0,5,0.1), target=1) {

  fit.rlm <- MASS::rlm(z~1, method="M", scale.est="MAD", psi=MASS::psi.huber)
  if (missing(sd0)) sd0 <- fit.rlm$s
  if (!is.numeric(sd0)) stop("Invalid sd0!")

  var0.gamma <- (sd0^2-1) / median(n.eff) / target
  c.val <- qnorm((1-level)/2, lower.tail=F) # critical value
  sd.supp <- sqrt(1 + var0.gamma*n.eff*target)
  lower0 <- -c.val*sd.supp # lower limits
  upper0 <- c.val*sd.supp # upper limits
  idx.null <- which(z>=lower0 & z<=upper0)
  z.null <- z[idx.null]
  N.null <- length(z.null)

  negloglkd <- function(pi.null, var.gamma) {
    sds <- sqrt(1 + var.gamma*n.eff*target)
    Q.nonnull <- pnorm(upper0[-idx.null]/sds[-idx.null]) -
      pnorm(lower0[-idx.null]/sds[-idx.null])
    negloglikelihood <-
      sum(log(sds[idx.null]) + z.null^2/(2*(sds[idx.null])^2)) -
      N.null * log(pi.null) - sum(log(1-pi.null*Q.nonnull))
    return(ifelse(is.nan(negloglikelihood), Inf, negloglikelihood))
  }

  mat <- sapply(grid.pi.null, function(pi.null) {
    mat <- sapply(grid.var.gamma, function(var.gamma) {
      c(pi.null, var.gamma,
        negloglkd(pi.null,var.gamma))
    })
    mat[,which.min(tail(mat, n=1))]
  })
  opt.init <- mat[,which.min(tail(mat, n=1))]
  inc.var.gamma <- min(diff(grid.var.gamma)) # incre of var.gamma
  grid.var.gamma <-
    seq(max(opt.init[2]-inc.var.gamma, 0),
        opt.init[2]+inc.var.gamma, inc.var.gamma/10)

  mat <- sapply(grid.var.gamma, function(var.gamma) {
    c(opt.init[1], var.gamma, negloglkd(opt.init[1], var.gamma))
  })
  opt.final <- mat[,which.min(tail(mat, n=1))]
  names(opt.final) <- c("pi.null", "var.gamma", "negloglkd.opt")
  return(opt.final)
}

#' Prepare Super Learner quantities
#'
#' @param object A fitted \code{drfDLNMadditive_fit} object with NB outcomes.
#' @param dat.c Required counterfactual data, including full exposure histories,
#'   on the same predictor scale expected by \code{predict()}. It must retain
#'   the fitted reference records and model structure.
#' @param verbose Whether to print progress.
#' @param AIC.only Used by \code{mDLNM::AICdata()} to prepare AIC quantities
#'   and a coefficient-sample variance instead of pairwise predictions.
#' @param nthreads Number of threads for deletion calculations, default \code{1}.
#'   Uses OpenMP for parallel calculation.
#' @param ... Unused.
#' @exportS3Method mDLNM::SLdata
SLdata.drfDLNMadditive_fit <- function(object, dat.c, verbose = FALSE,
                                     AIC.only = FALSE, nthreads = 1, ...) {
  if(missing(dat.c) || is.null(dat.c)) stop("dat.c is required.")
  if(length(list(...))) stop("Unused arguments.")
  if(identical(object$family, "gaussian")) stop("SLdata() requires negative binomial fits.")

  dat <- object$data
  point <- object$point
  SS.list <- dat$SS.list

  coef.samples <- NULL
  if(AIC.only) {
    samples <- object$CI.sample
    coef.samples <- cbind(samples$alpha_f_sample, samples$betaR_sample,
                          samples$betaF_sample)
    if(is.null(coef.samples) || nrow(coef.samples) < 2) {
      stop("AICdata() requires at least two stored coefficient draws in CI.sample.")
    }
  }

  mod.build <- drfDLNMadditivebuild(dat$y, dat$B_inner,
                    lapply(lapply(SS.list, "[[", "SSf"), "[[", "knots"),
                    lapply(lapply(SS.list, "[[", "SSw"), "[[", "knots"),
                    lapply(SS.list, "[[", "Sw"), lapply(SS.list, "[[", "SwR"),
                    lapply(SS.list, "[[", "Sf"), lapply(SS.list, "[[", "SfR"),
                    dat$Xrand, dat$Xfix, lapply(SS.list, "[[", "Zf"),
                    dat$offset$Xoffset, dat$r, point$alpha_f,
                    point$log_theta, point$log_smoothing_f, point$log_smoothing_w,
                    point$betaR, point$betaF, point$log_smoothing)

  mod.c <- stats::predict(object, dat = dat.c, .build.only = TRUE)
  out <- SLdrfDLNMadditive(mod.build$address.eigen, mod.c$address.eigen,
                          verbose, AIC.only, nthreads, coef.samples)
  theta <- exp(point$log_theta)
  if(AIC.only) {
    aic.names <- c("AIC_conventional", "AIC_cAIC")
    colnames(out$aic.loo) <- aic.names
    return(list(y = dat$y, row.id = rownames(object$modeldata), theta = theta,
                aic = stats::setNames(as.vector(out$aic), aic.names),
                aic.loo = out$aic.loo, delta = out$delta,
                delta.loo = as.vector(out$delta.loo), variance = stats::var(as.vector(out$delta.sample)),
                estimand = "ratio"))
  }
  return(list(y = dat$y, row.id = rownames(object$modeldata), theta = theta,
              mu.loo = as.vector(out$mu.loo), mu.lto = out$mu.lto,
              delta = out$delta, delta.loo = as.vector(out$delta.loo),
              estimand = "ratio"))
}

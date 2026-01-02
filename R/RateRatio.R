#' Compute Rate Ratio
#'
#' @param object An object.
#' @param ... Passed to methods.
#' @export
RateRatio <- function(object, ...) {
  UseMethod("RateRatio")
}


#' Title Rate Ratio for drfDLNMadditive_fit objects
#'
#' @param object object of class \code{drfDLNMadditive_fit}.
#' @param verbose whether to print messages during the process, default \code{FALSE}.
#'
#' @return
#' @importFrom mgcv s
#' @export
RateRatio.drfDLNMadditive_fit <- function(object, x0, x1, verbose = FALSE, ...) {


  pc <- object$data$pc
  kw <- object$data$kw
  kE <- object$data$kE
  sX.x.names.list <- object$data$sX.x.names.list
  B_inner <- object$data$B_inner
  M <- length(B_inner)
  maxL <- object$data$maxL

  R.CI <- nrow(object$CI.sample[[1]])
  n_alpha_f <- length(object$point$alpha_f)



  out.list <- lapply(1:M, function(i) {
    cen <- object$data$SS.list[[i]]$SSf$knots[as.integer(kE/2) + 2]
    if (verbose) cat("Exposure", i, "\n")
    ## point estimate
    alpha_f <- object$point$alpha_f[((i-1)*n_alpha_f/M + 1):(i*n_alpha_f/M)]
    Blag <- sapply(seq(0, maxL), function(l0) sidDLNM:::Bsplinevec2(l0, object$data$SS.list[[i]]$SSw$knots, 4))



    gridEl.x0 <- data.frame(E = x0[,i],
                            l = 0:maxL)
    gridEl.x1 <- data.frame(E = x1[,i],
                            l = 0:maxL)


    eta.E.0 <- sum(apply(gridEl.x0, 1, function(row.) {
      sidDLNM:::SurfaceEval(row.[1], cen, row.[2], alpha_f,
                            object$data$SS.list[[i]]$SSf$knots,
                            object$data$SS.list[[i]]$Zf,
                            Blag)
    }))

    eta.E.1 <- sum(apply(gridEl.x1, 1, function(row.) {
      sidDLNM:::SurfaceEval(row.[1], cen, row.[2], alpha_f,
                            object$data$SS.list[[i]]$SSf$knots,
                            object$data$SS.list[[i]]$Zf,
                            Blag)
    }))


    RR.est <- data.frame(eta.E.0 = eta.E.0,
                         eta.E.1 = eta.E.1,
                         RR = exp(eta.E.1 - eta.E.0))



    ## CI
    RR.sample <- lapply(1:R.CI, function(j) {
      if(verbose) {
        if(j %% 100 == 0) {
          cat("Sampling iteration:", j, " / ", R.CI, " for exposure",  i, "\n")
        }
      }
      alpha_f <- object$CI.sample$alpha_f_sample[j, ((i-1)*n_alpha_f/M + 1):(i*n_alpha_f/M)]

      eta.E.0 <- sum(apply(gridEl.x0, 1, function(row.) {
        sidDLNM:::SurfaceEval(row.[1], cen, row.[2], alpha_f,
                              object$data$SS.list[[i]]$SSf$knots,
                              object$data$SS.list[[i]]$Zf,
                              Blag)
      }))

      eta.E.1 <- sum(apply(gridEl.x1, 1, function(row.) {
        sidDLNM:::SurfaceEval(row.[1], cen, row.[2], alpha_f,
                              object$data$SS.list[[i]]$SSf$knots,
                              object$data$SS.list[[i]]$Zf,
                              Blag)
      }))


      return(data.frame(eta.E.0 = eta.E.0,
                           eta.E.1 = eta.E.1,
                           RR = exp(eta.E.1 - eta.E.0)))


    })


    RR.sample <- data.table::rbindlist(RR.sample)
    return(list(RR.est = RR.est,
                RR.sample = RR.sample))
  })

  out.list
}

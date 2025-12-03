#' Title fit a DRFDLNM using mgcv
#'
#' @param formula.list
#' @param kE
#' @param sXdat
#' @param M
#' @param verbose
#'
#' @returns
#' @importFrom tsModel Lag
#' @import mgcv
#' @import dlnm
#' @export
mgcv_fit <- function(formula.list, x.names.list, sXdat, kE, kw, verbose) {

  if(verbose) cat("Fitting a DRF-DLNM using mgcv package to obtain init ... \n")
  if(length(x.names.list) != 3) {
    error("not supported yet. TODO")
  }
  colnames(sXdat)[which(colnames(sXdat) == unlist(x.names.list))] <- c("x1", "x2", "x3")

  Q1 <- Lag(sXdat$x1, 0:maxL)
  Q2 <- Lag(sXdat$x2, 0:maxL)
  Q3 <- Lag(sXdat$x3, 0:maxL)


  cb1 <- crossbasis(Q1[,1],lag=c(0,maxL),argvar=list(fun='ps',df=kE-1),
                   arglag=list(fun='ps', df=kw))
  cb2 <- crossbasis(Q2[,1],lag=c(0,maxL),argvar=list(fun='ps',df=kE-1),
                    arglag=list(fun='ps', df=kw))
  cb3 <- crossbasis(Q3[,1],lag=c(0,maxL),argvar=list(fun='ps',df=kE-1),
                    arglag=list(fun='ps', df=kw))

  cbPengam1 <- cbPen(cb1)
  cbPengam2 <- cbPen(cb2)
  cbPengam3 <- cbPen(cb3)


  formula.list.mgcv <- Filter(Negate(is.null), formula.list[names(formula.list) %in% c("smooth","fe.varying")])


  formula.mgcv.pen <- "y~cb1+cb2+cb3"
  if(length(formula.list.mgcv) > 0) {
    formula.other.mgcv <- paste(formula.list.mgcv, collapse = "+")
    formula.other.mgcv <- gsub("~", "", formula.other.mgcv)
    formula.other.mgcv <- gsub("\"", "'", formula.other.mgcv)
    formula.mgcv.pen <- paste0(formula.mgcv.pen, "+", formula.other.mgcv)
  }
  formula.mgcv.pen <- as.formula(formula.mgcv.pen)

  # convert formula.mgcv.pen to string for printing
  formula.mgcv.pen.str <- deparse(formula.mgcv.pen)
  formula.mgcv.pen.str <- paste(formula.mgcv.pen.str, collapse = " ")
  if(verbose) cat("using mgcv::bam to fit model: ", formula.mgcv.pen.str, ". It might be slow, since some bam features are not supported for DRF-DLNMs.\n")

  # make bam function verbose
  b.opt <- bam(formula.mgcv.pen,
               family=nb(),
               data = sXdat,
               paraPen=list(cb1=cbPengam1,
                            cb2=cbPengam2,
                            cb3=cbPengam3),
               nthreads = 36,
               method = "REML")



  if(verbose) cat("DRF-DLNM using mgcv::bam(): Done! \n")

  return(list(mod.mgcv = b.opt)
  )
}






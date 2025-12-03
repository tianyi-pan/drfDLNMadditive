#' Title formula for Effect of exposure process
#'
#' @param x exposure variable
#' @param t time variable
#' @param by
#'
#' @return a formula for model fitting in \code{drfDLNMadditive}.
#' @export
sX <- function(t, x, by){
  xvar.add <- as.character(c(substitute(x)))
  xvar.list <- as.list(strsplit(xvar.add, split = "\\+")[[1]])
  xvar.list <- lapply(xvar.list, function(xx) gsub(" ", "", xx, fixed = TRUE))
  tvar <- as.character(substitute(t))
  if(missingArg(by)) {
    byvar <- NULL
  } else{
    byvar <- as.character(substitute(by))
  }
  return(list(x = xvar.list, t = tvar, by = byvar))
}

# function from mam https://github.com/awstringer1/mam/blob/master/R/helper.R
newsparsemat <- function(n,m) {
  methods::new('dgCMatrix',Dim=as.integer(c(n,m)),p=rep(0L,m+1),i=integer(0),x=numeric(0))
}


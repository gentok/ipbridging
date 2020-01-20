#' Estimating Normal Vector through Kernel Regularized Least Squares
#' 
#' @param xmat Matrix of OC coordinates (i.e., predictors).
#' @param resp Response Variable (i.e., ordered choices).
#' @param sample.size Sample size for data used for the estimation of KRLS.
#' @param ... Additional arguments passed to \code{\link[KRLS]{krls}}.
#' 
#' @return A vector of coefficients.
#' 
#' @seealso \code{\link[KRLS]{krls}}
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' This is a modified version of the code included in \code{\link[ooc]{ooc}} function.
#' 
#' @references 
#' Jeremy Ferwerda, Jens Hainmueller, Chad J. Hazlett (2017). Kernel-Based Regularized Least Squares in R (KRLS) and Stata (krls). Journal of Statistical Software, 79(3), 1-26. doi:10.18637/jss.v079.i03
#' 
#' @export

nv.krls <- function(xmat, resp, 
                    sample.size=300, ...) {
  
  temp.complete <- complete.cases(xmat & resp)
  X <- xmat[temp.complete,]
  y <- resp[temp.complete]
  select <- sample(1:length(y), sample.size, replace=TRUE)
  res <- KRLS::krls(X=X[select,], y=y[select], print.level=0, ...)
  coefs <- res$avgderivatives
  
  return(coefs)

}

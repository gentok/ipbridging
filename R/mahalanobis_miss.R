#' Mahalanobis distance (MD) for data with missing values
#'
#' This function is a modified version of \code{MDmiss} function in \code{modi} package.
#' 
#' For each observation the missing dimensions are omitted
#' before calculating the MD. The MD contains a correction
#' factor \eqn{p/q} to account for the number of observed values,
#' where \eqn{p} is the number of variables and \eqn{q} is the number of
#' observed dimensions for the particular observation.
#'
#' The function loops over the observations unless all values are non-missing.
#' This is not optimal if only a few missingness patterns occur. 
#' If no missing values occur
#' the function returns the Mahalanobis distance.
#'
#' @param data the data as a dataframe or matrix.
#' @param center the center to be used (may not contain missing values).
#' @param cov the covariance to be used (may not contain missing values).
#' @param ... Additional arguments passed to \link[stats]{mahalanobis} if no 
#' missing values exist.
#' @return The function returns a vector of the (squared) Mahalanobis distances.
#' @author Beat Hulliger (Original in \code{modi} package). Modified by Gento Kato.
#' @references Béguin, C., and Hulliger, B. (2004). Multivariate outlier detection
#' in incomplete survey data: The epidemic algorithm and transformed rank correlations.
#' Journal of the Royal Statistical Society, A167 (Part 2.), pp. 275-294.
#' @seealso \link[stats]{mahalanobis}
#' @export
mahalanobis_miss <- function(data, center, cov, ...) {
  
  # convert data frame to matrix if necessary
  if(!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  # get dimensions of matrix
  n <- nrow(data)
  p <- ncol(data)
  
  # returns TRUE if element of matrix is missing
  missings <- is.na(data)
  
  # count number of responded variables per observation
  resp.dim <- apply(!missings, 1, sum)
  
  # if no missing values, use stats::mahalanobis
  if (all(resp.dim)==p) {
    mdm <- stats::mahalanobis(data, center, cov, ...)
    return(mdm)
  }
  
  # center the data
  data <- sweep(data, 2, center)
  
  # create a numeric vector of length n to store results
  mdm <- numeric(n)
  
  # loop over observations
  for (i in 1:n) {
    
    # for obs. i, if there are any responses, then compute mahalanobis dist
    if (resp.dim[i] > 0) {
      
      resp <- !missings[i, ]
      x <- data[i, resp]
      mdm[i] <- t(x) %*% solve(cov[resp, resp]) %*% x
      
      # if no responses, return NA
    } else {
      
      mdm[i] <- NA
      
    }
  }
  
  # correction for number of responding variables
  mdm <- mdm * p / resp.dim
  
  return(mdm)
}
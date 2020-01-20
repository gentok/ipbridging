#' Bridging Two Ideal Point Estimates with OLS Mapping
#' 
#' @param ip1 Matrix or data.frame of 'reference' ideal points 
#' (i.e., \code{ip2} ideal points will be transformed and mapped to \code{ip1} space).
#' Rows are resepondents and columns are ideal point dimensions. 
#' The current code only allows 2 dimensions.
#' @param ip2 Matrix or data.frame of ideal points to be transformed.
#' @param anchorrows.ip1 Vector of row number of anchoring respondents in \code{ip1}.
#' @param anchorrows.ip2 Vector of row number of anchoring respondents in \code{ip2}.
#' Must be the same length as \code{anchorrows.ip1}.
#' 
#' @return A list with the following elements along with specified argument values:
#' \itemize{
#'   \item \code{ip1}: Original \code{ip1} matrix.
#'   \item \code{ip2_trans}: Transformed \code{ip2} matrix mapped on \code{ip1} space.
#'   \item \code{ip2_orig}: Original \code{ip2} matrix
#' }
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' 
#' @references Shor, B.; Berry, C. & McCarty, N. (2010), 'A Bridge to Somewhere: Mapping State and Congressional Ideology on a Cross-institutional Common Space', Legislative Studies Quarterly 35(3), 417--448. 
#' 
#' @import stats  
#'    
#' @export

bridge.olsmap <- function(ip1,
                          ip2,
                          anchorrows.ip1,
                          anchorrows.ip2) {
  
  ## Check inputs
  if (is.vector(ip1)) {
    warning("ip1 is a vector. Converted to matrix with one column.")
    ip1 <- as.matrix(ip1, ncol=1)
  }
  if (is.vector(ip2)) {
    warning("ip2 is a vector. Converted to matrix with one column.")
    ip2 <- as.matrix(ip2, ncol=1)
  }
  if (!is.matrix(ip1)&!is.data.frame(ip1)) stop("ip1 must be matrix or data.frame!") 
  if (!is.matrix(ip2)&!is.data.frame(ip2)) stop("ip2 must be matrix or data.frame!")
  if (ncol(ip1)!=ncol(ip2)) stop("ip1 and ip2 must have the same dimensions!")
  if (length(anchorrows.ip1)!=length(anchorrows.ip2)) {
    stop("anchorrows.ip1 and anchorrows.ip2 must have the same length!")
  }
  
  ## Generating Anchors 
  ac1 <- ip1[anchorrows.ip1,,drop=FALSE]
  ac2 <- ip2[anchorrows.ip2,,drop=FALSE]
  colnames(ac2) <- paste0("cd", 1:ncol(ac2))
  
  # Regression Formula
  f <- as.formula(paste("k ~", paste(colnames(ac2), collapse="+")))
  
  ## Models
  acmod <- plyr::alply(ac1, 2, function(k) lm(f, data=as.data.frame(cbind(k, ac2))))
  
  ## Prediction
  ip2df <- as.data.frame(ip2)
  colnames(ip2df) <- colnames(ac2)
  ip2_trans_f <- sapply(acmod, function(k) predict(k, ip2df))
  colnames(ip2_trans_f) <- colnames(ip2)
  
  ## For anchors, just use their ip1 coordinates
  ip2_trans_f[anchorrows.ip2,] <- ip1[anchorrows.ip1,]

  ######################
  ## Compiling Output ##
  ######################
  out <- list(ip1=ip1, 
              ip2_trans = ip2_trans_f,
              ip2_orig = ip2,
              anchorrows.ip1 = anchorrows.ip1,
              anchorrows.ip2 = anchorrows.ip2)
  
}  
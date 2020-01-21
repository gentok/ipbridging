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
  
  ## Exclude Missing Cases from the Analysis, if there are any
  comprows.ip1 <- which(complete.cases(ip1))
  comprows.ip2 <- which(complete.cases(ip2))
  ip1.rowid <- seq(1,nrow(ip1))
  ip2.rowid <- seq(1,nrow(ip2))
  if (length(c(comprows.ip1,comprows.ip2))<nrow(ip1)+nrow(ip2)) {
    
    # Store Originals
    ip1o <- ip1
    ip2o <- ip2
    
    # Dropping NAs from ip
    ip1 <- ip1o[comprows.ip1,]
    ip2 <- ip2o[comprows.ip2,]
    ip1.rowid <- ip1.rowid[comprows.ip1]
    ip2.rowid <- ip2.rowid[comprows.ip2]
    
    # Dropping NAs from anchorrows
    acmisloc <- unique(which(!anchorrows.ip1 %in% comprows.ip1), 
                       which(!anchorrows.ip2 %in% comprows.ip2))
    if (length(acmisloc)>0) {
      warning("Ideal points for some anchors are missing thus they are removed from anchors.")
      anchorrows.ip1 <- anchorrows.ip1o[-acmisloc]
      anchorrows.ip2 <- anchorrows.ip2o[-acmisloc]
    }
    
  } else {
    
    ip1o <- NULL
    ip2o <- NULL
    
  }
  
  ## Generating Anchors 
  ac1 <- ip1[which(ip1.rowid %in% anchorrows.ip1),,drop=FALSE]
  ac2 <- ip2[which(ip2.rowid %in% anchorrows.ip2),,drop=FALSE]
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
  ip2_trans_f[which(ip2.rowid %in% anchorrows.ip2),] <- 
    ip1[which(ip1.rowid %in% anchorrows.ip1),]
  
  # Putting NA Values Back In!
  if (!is.null(ip1o)) {
    
    ip1 <- ip1o
    ip2 <- ip2o
    ip2o[comprows.ip2,] <- ip2_trans_f
    ip2_trans_f <- ip2o
    
  }
  
  ######################
  ## Compiling Output ##
  ######################
  out <- list(ip1=ip1, 
              ip2_trans = ip2_trans_f,
              ip2_orig = ip2,
              anchorrows.ip1 = anchorrows.ip1,
              anchorrows.ip2 = anchorrows.ip2)
  
}  
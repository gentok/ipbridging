#' Bridging Two Ideal Point Estimates with Procrustes Transformation Method
#' 
#' @param ip1 Matrix or data.frame of 'reference' ideal points 
#' (i.e., \code{ip2} ideal points will be transformed and mapped to \code{ip1} space).
#' Rows are resepondents and columns are ideal point dimensions. 
#' The current code only allows 2 dimensions.
#' @param ip2 Matrix or data.frame of ideal points to be transformed.
#' @param anchorrows.ip1 Vector of row number of anchoring respondents in \code{ip1}.
#' @param anchorrows.ip2 Vector of row number of anchoring respondents in \code{ip2}.
#' Must be the same length as \code{anchorrows.ip1}.
#' @param opt If \code{TRUE}, conduct optimization of transformation through RANSAC
#' (random sample consensus). The default is \code{FALSE}.
#' @param opt.iter.n Number of iteration in the optimization of transformation 
#' matrix.
#' @param opt.sample.n Size of anchoring respondents to be sub-sampled at each iteration of 
#' optimization.
#' @param opt.th.inline Upper bound to determine inline respondents at each iteration of 
#' optimization. A respondent is considered 'inline' if the distance between transformed
#' \code{ip2} and \code{ip1} goes below this threshold.
#' @param blend If \code{FALSE}, do not use blending procedure. The default is \code{TRUE}.
#' @param blend.th1 If \code{blend==TRUE}, first threshold used in 'blending' procedure. If minimum difference 
#' between transformed \code{ip2} and \code{ip1} goes below this threshold, the 
#' transformed \code{ip2} is replaced with the closest \code{ip1}. 
#' @param blend.th2 If \code{blend==TRUE}, second threshold used in 'blending' procedure. If minimum difference 
#' between transformed \code{ip2} and \code{ip1} goes below this threshold (but above 
#' \code{blend.th1}), the transformed \code{ip2} is replaced with the value that 
#' intrapolates \code{ip1} that falls within this threshold. 
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
#' @export

bridge.procrustes <- function(ip1,
                              ip2,
                              anchorrows.ip1,
                              anchorrows.ip2, 
                              opt = FALSE,
                              opt.iter.n = 10000,
                              opt.sample.n = 30,
                              opt.th.inline = 0.5,
                              blend = TRUE,
                              blend.th1 = 0.05, 
                              blend.th2 = 0.15) {
  
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
  if (is.data.frame(ip1)) ip1 <- as.matrix(ip1)
  if (is.data.frame(ip2)) ip2 <- as.matrix(ip2)
  
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
    acmisloc <- unique(c(which(!anchorrows.ip1 %in% comprows.ip1), 
                         which(!anchorrows.ip2 %in% comprows.ip2)))
    if (length(acmisloc)>0) {
      warning("Ideal points for some anchors are missing thus they are removed from anchors.")
      anchorrows.ip1 <- anchorrows.ip1[-acmisloc]
      anchorrows.ip2 <- anchorrows.ip2[-acmisloc]
    }
    
  } else {
    
    ip1o <- NULL
    ip2o <- NULL
    
  }
  
  ## Generating Anchors 
  ac1 <- ip1[which(ip1.rowid %in% anchorrows.ip1),,drop=FALSE]
  ac2 <- ip2[which(ip2.rowid %in% anchorrows.ip2),,drop=FALSE]
  
  ##########################################
  ## Find Optimal Transformation Matrices ##
  ##########################################
  
  ## For 2-D Coordinates
  if (opt==TRUE) {
    
    ## Initial Values in Optimization
    num_in <- 0
    f_pool <- 0
    f_inl <- 0
    ite <- 1
    
    while (ite <= opt.iter.n) {
      
      ## Sampled Anchors for the simulation
      pool <- sample.int(nrow(ac2), size=opt.sample.n, replace=FALSE)
      ac1_r <- ac1[pool,]
      ac2_r <- ac2[pool,]
    
      #####################################
      ## Transforming sampled respondents #
      #####################################
      p <- MCMCpack::procrustes(ac2_r, ac1_r, translation=TRUE)

      ##############################################################
      ## Transforming the rest of (supposed) matching respondents ##
      ##############################################################
      ac2_trans <- p$s*(ac2%*%p$R) + (matrix(rep(1,nrow(ac2)), ncol=1) %*% t(p$tt))

      ## Find the number of inline respondents (under specified threshold)
      diff <- ac2_trans - ac1
      d2 <- diff^2
      euc2 <- apply(d2, 1, sum)
      eucd <- sqrt(euc2)
      t_num_in <- length(which(eucd <= opt.th.inline))
      inline <- which(eucd <= opt.th.inline)
      
      ## Update Output if more inline respondents found than previous iterations
      if(t_num_in > num_in){
        num_in <- t_num_in
        f_pool <- pool
        f_inl <- inline
      }
      
      ## Iteration Count
      ite <- ite + 1
      
    }
    
    #######################################
    ## Re-estimating p with all inliners ##
    #######################################
    r_p <- MCMCpack::procrustes(ac2[f_inl,], ac1[f_inl,], translation=TRUE)

    #########################################
    ## Transforming the "real" respondents ##
    #########################################
    ip2_trans <- r_p$s*(ip2%*%r_p$R) + (matrix(rep(1,nrow(ip2)), ncol=1) %*% t(r_p$tt))
    
  } else if (opt==FALSE) {
    
    ###############################
    ## Estimating Transformation ##
    ###############################
    r_p <- MCMCpack::procrustes(ac2, ac1, translation=TRUE)
    
    #############################
    ## Transforming All Points ##
    #############################
    ip2_trans <- r_p$s*(ip2%*%r_p$R) + (matrix(rep(1,nrow(ip2)), ncol=1) %*% t(r_p$tt))
    
  } else {
    
    stop("invalid 'opt' value!")
    
  }
  
  
  ## Create transformed coordinates ##
  if (blend==FALSE) {
    
    ip2_trans_f <- ip2_trans
    
  } else if (blend==TRUE) {
    
    ##############
    ## blending ##
    ##############
    
    ip2_blend <- function(j) {
      n_euc <- sqrt(apply(t(apply(ip1, 1, function(k) ip2_trans[j,] - k))^2, 1, sum))
      if (min(n_euc) < blend.th1|length(which(n_euc < blend.th2))==1) {
        
        # Return the coordinate in ip1 closest to the transformed coordinate
        return(ip1[which.min(n_euc),])
        
      } else {
        
        # Intrapolation if threshold is not met
        neighb <- which(n_euc < blend.th2)
        if (length(neighb)==0) {
          warning("No ip1 ideal point fall within difference threshold! NA generated.")
          return(rep(NA, ncol(ip1)))
        } else {
          neighb_p <- ip1[neighb,]
          weight <- n_euc[neighb]/sum(n_euc[neighb])
          return(apply(weight*neighb_p, 2, sum)) 
        }
        
      } 
      
    }
    ip2_trans_f <- matrix(NA, nrow=nrow(ip2), ncol=ncol(ip2))
    ## For non-anchors, set blended IPs
    ip2_trans_f[-which(ip2.rowid %in% anchorrows.ip2),] <- 
      do.call("rbind", lapply(seq(1,nrow(ip2))[-which(ip2.rowid %in% anchorrows.ip2)], ip2_blend) )
    ## For anchors, just use their ip1 coordinates
    ip2_trans_f[which(ip2.rowid %in% anchorrows.ip2),] <- 
      ip1[which(ip1.rowid %in% anchorrows.ip1),]
    
  } else {
    
    stop("invalid 'blend' value!")
    
  }
  
  ################################
  ## Putting NA Values Back In! ##
  ################################
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
              anchorrows.ip2 = anchorrows.ip2, 
              opt = opt,
              opt.iter.n = opt.iter.n,
              opt.sample.n = opt.sample.n,
              opt.th.inline = opt.th.inline,
              blend.th1 = blend.th1, 
              blend.th2 = blend.th2)
  
}

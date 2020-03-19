#' Bridging Two Ideal Point Estimates with Linear Transformation Method
#' 
#' @param ip1 Matrix or data.frame of 'reference' ideal points 
#' (i.e., \code{ip2} ideal points will be transformed and mapped to \code{ip1} space).
#' Rows are resepondents and columns are ideal point dimensions. 
#' The current code only allows 2 dimensions.
#' @param ip2 Matrix or data.frame of ideal points to be transformed.
#' @param anchorrows.ip1 Vector of row number of anchoring respondents in \code{ip1}.
#' @param anchorrows.ip2 Vector of row number of anchoring respondents in \code{ip2}.
#' Must be the same length as \code{anchorrows.ip1}.
#' @param method Method of bridging. Currently, following methods are aviailable:
#' \itemize{
#'   \item \code{"procrustes"} (default): Procrustes transformation method. 
#'   Based on anchor cases, this method provides restricted non-parametric procedure to 
#'   find optimal transformation matrix to bridge ideal point estimates.
#'   \item \code{"homography"}: Homography transformation method. 
#'   Based on anchor cases, this method provides non-parametric procedure to 
#'   find optimal transformation matrix to bridge ideal point estimates.
#'   \item \code{"olsmap"}: OLS mapping method, Based on anchor cases, 
#'   use OLS regression to map \code{d2} ideal point coordinates on 
#'   \code{d1} ideal point space.  
#' }
#' @param trans.ip2 If \code{TRUE} (default), transform \code{ip2} to map them on \code{ip1} space.
#' If \code{FALSE}, transform \code{ip1} to map them on \code{ip2} space.
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
#'   \item \code{ip2_trans} or \code{ip1_trans}: Transformed matrix.
#'   \item \code{ip1}: Original \code{ip1} matrix.
#'   \item \code{ip2}: Original \code{ip2} matrix
#' }
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' 
#' @export

bridge.linearmap <- function(ip1,
                             ip2,
                             anchorrows.ip1,
                             anchorrows.ip2,
                             method = "procrustes",
                             trans.ip2 = TRUE,
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
  
  ## For now, the homography method is only applicable to 
  ## 2-D ideal point coordinates
  if (method=="homography") {
    if (ncol(ip1)!=2) stop("The homography method is currently only applicable to 2-D coordinates!")
  }
  
  ## If transforming ip1 instead of ip2...
  if (trans.ip2==FALSE) {
    
    # Replacing Ideal Points
    tmp <- ip1
    ip1 <- ip2
    ip2 <- tmp
    
    # Replacing anchorrows
    tmp <- anchorrows.ip1
    anchorrows.ip1 <- anchorrows.ip2
    anchorrows.ip2 <- tmp
    
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
      
      if (method=="procrustes") {
        
        ######################################
        ## Transforming sampled respondents ##
        ######################################
        p <- MCMCpack::procrustes(ac2_r, ac1_r, translation=TRUE)
        
        ##############################################################
        ## Transforming the rest of (supposed) matching respondents ##
        ##############################################################
        ac2_trans <- p$s*(ac2%*%p$R) + (matrix(rep(1,nrow(ac2)), ncol=1) %*% t(p$tt))

      } else if (method=="homography") {
        
        ######################################
        ## Transforming sampled respondents ##
        ######################################
        A <- do.call("rbind", 
                     lapply(1:nrow(ac2_r), 
                            function(i) matrix(c(ac2_r[i,],1,0,0,0,
                                                 -ac2_r[i,1]*ac1_r[i,1],
                                                 -ac2_r[i,2]*ac1_r[i,1],
                                                 0,0,0,ac2_r[i,],1,
                                                 -ac2_r[i,1]*ac1_r[i,2],
                                                 -ac2_r[i,2]*ac1_r[i,2]), 
                                               ncol=8, byrow = TRUE)
                     ))
        b <- unlist(lapply(1:nrow(ac2_r), function(i) ac1_r[i,]))
        h <- solve(t(A)%*%A)%*%t(A)%*%b
        #h <- matrix(c(h,1), nrow=3, byrow=TRUE)
        #h <- solve(h)
        #h <- matrix(t(h), ncol=1, byrow=FALSE)[1:8,]
        
        ##############################################################
        ## Transforming the rest of (supposed) matching respondents ##
        ##############################################################
        ac2_trans <- 
          do.call("rbind",
                  lapply(1:nrow(ac2),
                         function(i) 
                           matrix(c((h[1]*ac2[i,1] + h[2]*ac2[i,2] + h[3])/
                                      (h[7]*ac2[i,1] + h[8]*ac2[i,2] + 1), 
                                    (h[4]*ac2[i,1] + h[5]*ac2[i,2] + h[6])/
                                      (h[7]*ac2[i,1] + h[8]*ac2[i,2] + 1)),
                                  nrow=1)
                  ))

      } else if (method=="olsmap") {
        
        #####################################
        ## Transforming sampled respondents #
        #####################################
        # Regression Formula
        f <- as.formula(paste("k ~", paste(colnames(ac2_r), collapse="+")))
        ## Models
        acmod <- plyr::alply(ac1_r, 2, 
                             function(k) lm(f, data=as.data.frame(cbind(k, ac2_r))))
        
        ##############################################################
        ## Transforming the rest of (supposed) matching respondents ##
        ##############################################################
        ## Prediction
        ac2df <- as.data.frame(ac2)
        colnames(ac2df) <- colnames(ac2)
        ac2_trans <- sapply(acmod, function(k) predict(k, ac2df))
        colnames(ac2_trans) <- colnames(ac2)

      } else {
        
        stop("invalid 'method' value!")
      }

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
    
    if (method=="procrustes") {
      
      #######################################
      ## Re-estimating p with all inliners ##
      #######################################
      r_p <- MCMCpack::procrustes(ac2[f_inl,], ac1[f_inl,], translation=TRUE)
      
      #########################################
      ## Transforming the "real" respondents ##
      #########################################
      ip2_trans <- r_p$s*(ip2%*%r_p$R) + (matrix(rep(1,nrow(ip2)), ncol=1) %*% t(r_p$tt))
      
    } else if (method=="homography") {
      
      #######################################
      ## Re-estimating h with all inliners ##
      #######################################
      r_A <- do.call("rbind", 
                     lapply(1:num_in, 
                            function(i) matrix(c(ac2[f_inl[i],],1,0,0,0,
                                                 -ac2[f_inl[i],1]*ac1[f_inl[i],1],
                                                 -ac2[f_inl[i],2]*ac1[f_inl[i],1],
                                                 0,0,0,ac2[f_inl[i],],1,
                                                 -ac2[f_inl[i],1]*ac1[f_inl[i],2],
                                                 -ac2[f_inl[i],2]*ac1[f_inl[i],2]), 
                                               ncol=8, byrow=TRUE)
                     ))
      r_b <- unlist(lapply(1:num_in, function(i) ac1[f_inl[i],]))
      r_h <- solve(t(r_A)%*%r_A)%*%t(r_A)%*%r_b
      #r_h <- matrix(c(r_h,1), nrow=3, byrow=TRUE)
      #r_h <- solve(r_h)
      #r_h <- matrix(t(r_h), ncol=1, byrow=FALSE)[1:8,]
      
      #########################################
      ## Transforming the "real" respondents ##
      #########################################
      ip2_trans <- 
        do.call("rbind",
                lapply(1:nrow(ip2),
                       function(i) 
                         matrix(c((r_h[1]*ip2[i,1] + r_h[2]*ip2[i,2] + r_h[3])/
                                    (r_h[7]*ip2[i,1] + r_h[8]*ip2[i,2] +1), 
                                  (r_h[4]*ip2[i,1] + r_h[5]*ip2[i,2] + r_h[6])/
                                    (r_h[7]*ip2[i,1] + r_h[8]*ip2[i,2] + 1)),
                                nrow=1)
                ))
      
    } else if (method=="olsmap") {
      
      ###########################################
      ## Re-estimating model with all inliners ##
      ###########################################
      # Regression Formula
      r_f <- as.formula(paste("k ~", paste(colnames(ac2[f_inl,]), collapse="+")))
      ## Models
      r_acmod <- plyr::alply(ac1[f_inl,], 2, 
                             function(k) lm(r_f, data=as.data.frame(cbind(k, ac2[f_inl,]))))
      
      #########################################
      ## Transforming the "real" respondents ##
      #########################################
      ## Prediction
      ip2df <- as.data.frame(ip2)
      colnames(ip2df) <- colnames(ac2)
      ip2_trans <- sapply(r_acmod, function(k) predict(k, ip2df))
      colnames(ip2_trans) <- colnames(ip2)

    } else {
      
      stop("invalid 'method' value!")
    }

  } else if (opt==FALSE) {
    
    if (method=="procrustes") {
      
      ###############################
      ## Estimating Transformation ##
      ###############################
      r_p <- MCMCpack::procrustes(ac2, ac1, translation=TRUE)
      
      #############################
      ## Transforming All Points ##
      #############################
      ip2_trans <- r_p$s*(ip2%*%r_p$R) + (matrix(rep(1,nrow(ip2)), ncol=1) %*% t(r_p$tt))

    } else if (method=="homography") {
      
      ###################################
      ## Estimating h with all anchors ##
      ###################################
      r_A <- do.call("rbind", 
                     lapply(1:num_in, 
                            function(i) matrix(c(ac2[i,],1,0,0,0,
                                                 -ac2[i,1]*ac1[i,1],
                                                 -ac2[i,2]*ac1[i,1],
                                                 0,0,0,ac2[i,],1,
                                                 -ac2[i,1]*ac1[i,2],
                                                 -ac2[i,2]*ac1[i,2]), 
                                               ncol=8, byrow=TRUE)
                     ))
      r_b <- unlist(lapply(1:num_in, function(i) ac1[f_inl[i],]))
      r_h <- solve(t(r_A)%*%r_A)%*%t(r_A)%*%r_b
      #r_h <- matrix(c(r_h,1), nrow=3, byrow=TRUE)
      #r_h <- solve(r_h)
      #r_h <- matrix(t(r_h), ncol=1, byrow=FALSE)[1:8,]
      
      #########################################
      ## Transforming the "real" respondents ##
      #########################################
      ip2_trans <- 
        do.call("rbind",
                lapply(1:nrow(ip2),
                       function(i) 
                         matrix(c((r_h[1]*ip2[i,1] + r_h[2]*ip2[i,2] + r_h[3])/
                                    (r_h[7]*ip2[i,1] + r_h[8]*ip2[i,2] +1), 
                                  (r_h[4]*ip2[i,1] + r_h[5]*ip2[i,2] + r_h[6])/
                                    (r_h[7]*ip2[i,1] + r_h[8]*ip2[i,2] + 1)),
                                nrow=1)
                ))

    } else if (method=="olsmap") {
      
      ###############################
      ## Estimating Transformation ##
      ###############################
      # Regression Formula
      f <- as.formula(paste("k ~", paste(colnames(ac2), collapse="+")))
      ## Models
      acmod <- plyr::alply(ac1, 2, function(k) lm(f, data=as.data.frame(cbind(k, ac2))))
      
      #############################
      ## Transforming All Points ##
      #############################
      ## Prediction
      ip2df <- as.data.frame(ip2)
      colnames(ip2df) <- colnames(ac2)
      ip2_trans <- sapply(acmod, function(k) predict(k, ip2df))
      colnames(ip2_trans) <- colnames(ip2)

    } else {
      
      stop("invalid 'method' value!")
    }

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
    # Transformed Coordinates
    ip2_trans_f <- matrix(NA, nrow=nrow(ip2), ncol=ncol(ip2))
    ip2_trans_f <- do.call("rbind", lapply(seq(1,nrow(ip2)), ip2_blend))
    # ## For non-anchors, set blended IPs
    # ip2_trans_f[-which(ip2.rowid %in% anchorrows.ip2),] <- 
    #   do.call("rbind", lapply(seq(1,nrow(ip2))[-which(ip2.rowid %in% anchorrows.ip2)], ip2_blend) )
    # ## For anchors, just use their ip1 coordinates
    # ip2_trans_f[which(ip2.rowid %in% anchorrows.ip2),] <- 
    #   ip1[which(ip1.rowid %in% anchorrows.ip1),]
    
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
  if (trans.ip2) {
    
    out <- list(ip2_trans = ip2_trans_f,
                ip1=ip1, 
                ip2=ip2,
                anchorrows.ip1 = anchorrows.ip1,
                anchorrows.ip2 = anchorrows.ip2,
                method = method,
                opt = opt,
                opt.iter.n = opt.iter.n,
                opt.sample.n = opt.sample.n,
                opt.th.inline = opt.th.inline,
                blend.th1 = blend.th1, 
                blend.th2 = blend.th2)

  } else {
    
    out <- list(ip1_trans = ip2_trans_f,
                ip1=ip2, 
                ip2=ip1,
                anchorrows.ip1 = anchorrows.ip2,
                anchorrows.ip2 = anchorrows.ip1,
                method = method,
                opt = opt,
                opt.iter.n = opt.iter.n,
                opt.sample.n = opt.sample.n,
                opt.th.inline = opt.th.inline,
                blend.th1 = blend.th1, 
                blend.th2 = blend.th2)
    
  }
  
  return(out)
  
}

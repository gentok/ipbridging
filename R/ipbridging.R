#' Bridging Ideal Point Estimates from Two Data Sets
#' 
#' @description Toolbox function that bridges ideal point estimates computed from 
#' two separate datasets. 
#' 
#' @param d1 Matrix or data.frame of 'reference' responses.
#' (i.e., In mapping methods, ideal points based on \code{d2} will be 
#' transformed and mapped to ideal point space based on \code{d1}).
#' Rows are resepondents and columns are questions/issues. 
#' Alternatively, if \code{input.type=="idealpoints"}, one can directly 
#' set matrix or data.frame of ideal points (rows are respondents and 
#' columns are ideal point dimensions).
#' @param d2 Matrix or data.frame of responses (or ideal points) to be transformed.
#' @param bridge.method Method of bridging. Currently, following methods are aviailable:
#' \itemize{
#'   \item \code{"procrustes"} (default): Procrustes transformation method. 
#'   Based on anchor cases, this method provides restricted non-parametric procedure to 
#'   find optimal transformation matrix to bridge ideal point estimates (calls \code{\link{bridge.homography}}).
#'   \item \code{"homography"}: Homography transformation method. 
#'   Based on anchor cases, this method provides non-parametric procedure to 
#'   find optimal transformation matrix to bridge ideal point estimates (calls \code{\link{bridge.homography}}).
#'   \item \code{"joint"}: Joint scaling method. Simple method to combine \code{d1} and 
#'   \code{d2} and jointly compute ideal points. The procedure described in Jessee (2016).
#'   \item \code{"modelmap"}: Use the model estimated on \code{d1} to map
#'   ideal point coordinates of \code{d2} respondents. This method does not require 
#'   anchor cases. Can be used only if \code{input.type=="responses"} and 
#'   \code{ip.method=="irtMCMC"}. See Jessee (2016).
#'   \item \code{"olsmap"}: OLS mapping method, Based on anchor cases, 
#'   use OLs regression to map \code{d2} ideal point coordinates on 
#'   \code{d1} ideal point space. The procedure described in Shor et al. (2010). 
#'   \item \code{"anchoredpooling"}: Use anchor cases to pool \code{d1} and \code{d2} 
#'   and jointly compute ideal points. Issues/questions in two datasets are assumed 
#'   to be independent from each other. See Shor et al. (2010) for more details.
#' }
#' @param anchors.method Method to select anchors to bridge ideal points. 
#' Following values are currently available. 
#' Ignored if \code{bridge.method} is \code{"joint"} or \code{"modelmap"}.
#' \itemize{
#'   \item \code{"subsample"} (default): Sample subset of \code{d1} and \code{d2} and 
#'   use them as anchors. Sampled anchors from another dataset are added to the 
#'   target dataset when computing ideal points on each dataset. 
#'   See \code{anchors.subsample.pr} and \code{anchors.subsample.method} for 
#'   sampling options.
#'   \item \code{"selectrows"}: Use this option if specific anchors are already 
#'   included in both \code{d1} and \code{d2}. Row numbers for anchors must be set 
#'   in \code{anchors.selectrows.d1} and \code{anchors.selectrows.d2}, respectively
#'   (Two vectors must have the same length and each element must corresponds to 
#'   the same anchor respondent). This is the only option for anchors 
#'   selection if \code{input.type=="idealpoints"}. 
#'   \item \code{"newdata"} Use new dataset that only contains anchors. The data.frame
#'   assigned in \code{anchors.newdata} are appended to both \code{d1} and \code{d2}. 
#' }
#' Note that \code{"subsample"} and \code{"newdata"} requires \code{d1} and \code{d2}
#' to have identical columns. This is not required for \code{"selectrows"}. 
#' @param anchors.subsample.method Method to sample anchors. Following values 
#' are currently available: 
#' \itemize{
#'   \item \code{"random"} (default): randomly sample subset of both datasets and 
#'   use them as anchors. The sample size is determined by \code{anchors.subsample.pr}.
#'   \item \code{"random.d1"}: randomly sample subset of \code{d1} and 
#'   use them as anchors. The sample size is determined by \code{anchors.subsample.pr}.
#'   \item \code{"random.d2"}: randomly sample subset of \code{d2} and 
#'   use them as anchors. The sample size is determined by \code{anchors.subsample.pr}.
#'   \item \code{"selectrows"}: Use this option if specific anchors are already known (
#'   but not shared between datasets). Row numbers for anchors must be set 
#'   in \code{anchors.selectrows.d1} and \code{anchors.selectrows.d2}, respectively
#'   (Two vectors does not have to have the same length).
#' }
#' @param anchors.subsample.pr Used when  \code{anchors.method=="subsample"}. 
#' Proportion of dataset (i.e., \code{d1} and \code{d2}) 
#' to be sampled as anchor. The propotion is determined by the dataset of smaller size if 
#' \code{anchors.subsample.method=="random"}, specific dataset's size if 
#' \code{anchors.subsample.method} is \code{"random.d1"} or \code{"random.d2"}. 
#' @param anchors.subsample.wgt.d1 Optional weights passed to \code{prob} option of 
#' \code{\link[base]{sample}} function used in random sampling of anchors from \code{d1} when 
#' \code{anchors.subsample.method} is \code{"random"} or \code{"random.d1"}. If not \code{NULL}, 
#' it must be a vector of length equals to \code{nrow(d1)}. 
#' @param anchors.subsample.wgt.d2 Optional weights passed to \code{prob} option of 
#' \code{\link[base]{sample}} function used in random sampling of anchors from \code{d2} when 
#' \code{anchors.subsample.method} is \code{"random"} or \code{"random.d2"}. If not \code{NULL}, 
#' it must be a vector of length equals to \code{nrow(d1)}. 
#' @param anchors.selectrows.d1 Must be the vector of positive integers.  
#' Specify the row numbers of anchors in \code{d1}. 
#' If \code{anchors.method=="selectrows"}, must be the same length and with perfect correspondence with \code{anchors.selectrows.d2}.
#' If \code{anchors.method=="subsample"} and \code{anchors.subsample.method=="selectrows"}, 
#' the above condition is not required. 
#' @param anchors.selectrows.d2 Must be the vector of positive integers.  
#' Specify the row numbers of anchors in \code{d2}.
#' @param anchors.newdata Matrix or data.frame of anchors. Must have identical 
#' columns as \code{d1} and \code{d2}. Used if \code{anchors.method=="newdata"}. 
#' @param input.type The input type of \code{d1} and \code{d2}. If \code{"responses"} (default), 
#' dataset values are responses to multiple issues or votes. 
#' If \code{"idealpoints"}, dataset values are ideal points.   
#' @param ip.method Method of ideal point computation (used only when \code{input.type=="response"}. Following values are currently available.
#' \itemize{
#'   \item \code{"ooc"} (default): Ordered optimal classification. Calls \code{\link{oocflex}} function.
#'   \item \code{"blackbox"}: Blackbox scaling. Calls \code{\link[basicspace]{blackbox}} function from \code{blackbox} package (\code{minscale=10} by default).
#'   \item \code{"oc"}: Optimal classification. Calls \code{\link[oc]{oc}} function from \code{oc} package (\code{polarity=c(1,1)} by default). 
#'   \item \code{"wnominate"}: W-NOMINATE. Calls \code{\link[wnominate]{wnominate}} function from \code{wnominate} package  (\code{polarity=c(1,1)} by default).
#'   \item \code{"irtMCMC"}: IRT model via Markov chain Monte Carlo method. Calls \code{\link[pscl]{ideal}} function from \code{pscl} package.
#' }
#' For \code{ip.method=="ooc"}, the inputs \code{d1} and \code{d2} must be positive integer. 
#' For all other methods, the inputs must be a dummy variable indication 0=Nay and 1=Yea (and 9 for unit nonresponse).
#' All methods allow item non-response as NA.
#' @param ip.dims Number of dimension in ideal point computation. Must be a positive 
#' integer between 1 and 10. If \code{bridge.method=="homography"}, only 2 is accepted 
#' at this point.  
#' @param ip.polarity.d1 A vector specifying the row number of the \code{d1}(or pooled data)'s respondent(s) 
#' constrained to have a positive (i.e., right-wing or conservative) score on each dimension.
#' Used when \code{ip.method} is \code{"ooc"}, \code{"oc"}, or \code{"wnominate"}. 
#' @param ip.polarity.d2 A vector specifying the row number of the \code{d2}'s respondent(s) 
#' constrained to have a positive (i.e., right-wing or conservative) score on each dimension. 
#' Used when \code{ip.method} is \code{"ooc"}, \code{"oc"}, or \code{"wnominate"}. 
#' @param tr.opt Used if \code{bridge.method} is a transformation method (i.e., \code{"procrustes"}, \code{"homography"}, or \code{"olsmap"}). 
#' If \code{TRUE}, conduct optimization of transformation through RANSAC
#' (random sample consensus). The default is \code{FALSE}.
#' @param tr.opt.iter.n Used if \code{bridge.method} is a transformation method (i.e., \code{"procrustes"}, \code{"homography"}, or \code{"olsmap"}). 
#' Number of iteration in the optimization of transformation matrix.
#' @param tr.opt.sample.n Used if \code{bridge.method} is a transformation method (i.e., \code{"procrustes"}, \code{"homography"}, or \code{"olsmap"}). 
#' Size of anchoring respondents to be sub-sampled at each iteration of 
#' optimization.
#' @param tr.opt.th.inline Used if \code{bridge.method} is a transformation method (i.e., \code{"procrustes"}, \code{"homography"}, or \code{"olsmap"}). 
#' Upper bound to determine inline respondents at each iteration of 
#' optimization. A respondent is considered 'inline' if the distance between transformed
#' \code{ip2} and \code{ip1} goes below this threshold.
#' @param tr.blend Used if \code{bridge.method} is a transformation method (i.e., \code{"procrustes"}, \code{"homography"}, or \code{"olsmap"}). 
#' If \code{FALSE}, do not use blending procedure in homograhy transformation. 
#' The default is \code{TRUE}.
#' @param tr.blend.th1 Used if \code{bridge.method} is a transformation method (i.e., \code{"procrustes"}, \code{"homography"}, or \code{"olsmap"}). 
#' First threshold used in 'blending' procedure. If minimum difference 
#' between transformed \code{ip2} and \code{ip1} goes below this threshold, the 
#' transformed \code{ip2} is replaced with the closest \code{ip1}. 
#' @param tr.blend.th2 Used if \code{bridge.method} is a transformation method (i.e., \code{"procrustes"}, \code{"homography"}, or \code{"olsmap"}). Second threshold used in 'blending' procedure. If minimum difference 
#' between transformed \code{ip2} and \code{ip1} goes below this threshold (but above 
#' \code{blend.th1}), the transformed \code{ip2} is replaced with the value that 
#' intrapolates \code{ip1} that falls within this threshold. 
#' @param ... Additional arguments passed to the ideal point estimation function 
#' called in \code{ip.method}.
#' 
#' @return A list with the following elements along with specified argument values
#' \itemize{
#'   \item \code{ip1}: Original \code{ip1} matrix.
#'   \item \code{ip2_trans}: Transformed \code{ip2} matrix mapped on \code{ip1} space.
#'   \item \code{ip2_orig}: Original \code{ip2} matrix
#' }
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' 
#' @seealso \code{\link{ipest}}, \code{\link{setanchors}}, 
#' \code{\link{bridge.procrustes}}, \code{\link{bridge.homography}}, 
#' \code{\link{bridge.olsmap}},
#' \code{\link{oocflex}}, \code{\link[oc]{oc}}, \code{\link[basicspace]{blackbox}},
#' \code{\link[wnominate]{wnominate}}, \code{\link[pscl]{ideal}}
#' 
#' @references 
#' \itemize{
#'   \item Jessee, S. (2016), '(How) Can We Estimate the Ideology of Citizens and Political Elites on the Same Scale?', American Journal of Political Science 60(4), 1108--1124.
#'   \item Shor, B.; Berry, C. & McCarty, N. (2010), 'A Bridge to Somewhere: Mapping State and Congressional Ideology on a Cross-institutional Common Space', Legislative Studies Quarterly 35(3), 417--448.
#'   
#' }
#' 
#' @import stats
#' @import oc
#' @import ooc
#' @import wnominate
#' @import pscl
#' @import basicspace
#' @import emIRT
#' 
#' @export

ipbridging <- function(d1, d2, 
                       bridge.method = "procrustes", 
                       anchors.method = "subsample",
                       anchors.subsample.method = "random",
                       anchors.subsample.pr = 0.1,
                       anchors.subsample.wgt.d1 = NULL,
                       anchors.subsample.wgt.d2 = NULL,
                       anchors.selectrows.d1 = 1:100,
                       anchors.selectrows.d2 = 1:100,
                       anchors.newdata = NULL,
                       input.type="responses",
                       ip.method="ooc", 
                       ip.dims=2,
                       ip.polarity.d1 = c(1,1),
                       ip.polarity.d2 = c(1,1),
                       tr.opt = FALSE,
                       tr.opt.iter.n = 10000,
                       tr.opt.sample.n = 30,
                       tr.opt.th.inline = 0.5,
                       tr.blend = TRUE,
                       tr.blend.th1 = 0.05, 
                       tr.blend.th2 = 0.15,
                       ...) {
  
  ## Check Inputs
  if (is.vector(d1)) {
    warning("d1 is a vector. Converted to matrix with one column.")
    d1 <- as.matrix(d1, ncol=1)
  }
  if (is.vector(d2)) {
    warning("d2 is a vector. Converted to matrix with one column.")
    d2 <- as.matrix(d2, ncol=1)
  }
  if (!is.matrix(d1)&!is.data.frame(d1)) stop("d1 must be matrix or data.frame!") 
  if (!is.matrix(d2)&!is.data.frame(d2)) stop("d2 must be matrix or data.frame!") 
  if (is.data.frame(d1)) d1 <- as.matrix(d1)
  if (is.data.frame(d2)) d2 <- as.matrix(d2)
  
  # Respondents Data
  respdt <- data.frame(allid = seq(1,nrow(d1)+nrow(d2)),
                       subid = c(seq(1,nrow(d1)),seq(1,nrow(d2))),
                       data = c(rep(1, nrow(d1)), rep(2, nrow(d2))),
                       isanchor = 0)
  
  cat(paste0("\nBridging ideal points through ", bridge.method," method:\n\n"))
  
  ## Joint scaling As a Bridging Method
  if (bridge.method=="joint") {

    if (ncol(d1)!=ncol(d2)) stop("d1 and d2 must have the same number of columns!")
    
    if (input.type=="idealpoints") {
      
      stop('input.type=="idealpoints" is incompatible with bridge.method=="joint"!')
      
    } else if (input.type=="responses") {
      
      cat("Generating bridged ideal points...\n")
      
      # Anchor Identifier (NA in joint)
      respdt$isanchor <- NA
      
      # Estimate Ideal Points
      ip_joint_est <- ipest(rbind(d1,d2), method=ip.method, 
                            dims=ip.dims, polarity=ip.polarity.d1, ...)
      ip_joint_f <- ip_joint_est$ip_f
      ip_joint <- ip_joint_est$ip
      
    } else {
      
      stop("Invalid 'input.type' value!")
      
    }
    
    ## Add IP values to respdt
    for(i in 1:ncol(ip_joint)) respdt[,paste0("bridged",i,"D")] <- ip_joint[,i]
    for(i in 1:ncol(ip_joint)) respdt[,paste0("ip1_coord",i,"D")] <- NA
    for(i in 1:ncol(ip_joint)) respdt[,paste0("ip2_coord",i,"D")] <- NA
    
    cat("DONE!\n\n")
    
    ## Compile Output
    out <- list(bridge.data = respdt, 
                bridge.model = ip_joint_f,
                bridge.method = bridge.method, 
                ip.model.d1 = NULL,
                ip.model.d2 = NULL,
                ip.method = ip.method,
                ip.dims = ip.dims,
                input.type = input.type,
                anchors.method = NULL)
    return(out)
    
  ## All Other Bridging Methods using Linear Mapping with Anchors    
  } else if (bridge.method%in%c("procrustes","homography","olsmap")) {
    
    ## Generate Ideal Point Estimates from Two Separate Data Sets
    if (input.type=="idealpoints") {
      
      if (ncol(d1)!=ncol(d2)) stop("d1 and d2 must have the same number of columns!")
      
      # Input is already ideal points
      ip1 <- d1
      ip2 <- d2
      
      # Update ip relevant parameters
      ip.method <- ip1_f <- ip2_f <- NULL
      ip.dims <- ncol(ip1)
      
      # Setting Anchors
      cat("Setting anchors...\n\n")
      if (anchors.method!="selectrows") {
        
        warning("anchors.method must be 'selectrows' for input.type=='idealpoints'. Changed anchors.method to 'selectrows'.")
        anchors.method = "selectrows"
        
      }
      if (anchors.method=="selectrows") {
        
        if (length(anchors.selectrows.d1)!=length(anchors.selectrows.d2)) {
          stop("anchors.selectrows.d1 and anchors.selectrows.d2 must have the same length!")
        }
        
        # Anchor Ideal Points
        anchorrows.ip1 <- anchors.selectrows.d1
        anchorrows.ip2 <- anchors.selectrows.d2

      }
      
      
    } else if (input.type=="responses") {
      
      outac <- setanchors(d1, d2,
                          anchors.method,# = "subsample",
                          anchors.subsample.method,# = "random",
                          anchors.subsample.pr,# = 0.1,
                          anchors.subsample.wgt.d1,# = NULL,
                          anchors.subsample.wgt.d2,# = NULL,
                          anchors.selectrows.d1,# = 1:100,
                          anchors.selectrows.d2,# = 1:100,
                          anchors.newdata)# = NULL)
      d1x <- outac$d1x
      d2x <- outac$d2x
      respdt <- outac$respdt
      anchorrows.ip1 <- outac$anchorrows.ip1
      anchorrows.ip2 <- outac$anchorrows.ip2
      
      ## Compute Ideal Points
      cat("Generating ideal points on subsets...\n")
      # Ideal Point 1
      ip1_est <- ipest(d1x, method=ip.method, 
                            dims=ip.dims, polarity=ip.polarity.d1, ...)
      ip1_f <- ip1_est$ip_f
      ip1 <- ip1_est$ip
      # Ideal Point 2
      ip2_est <- ipest(d2x, method=ip.method, 
                       dims=ip.dims, polarity=ip.polarity.d2, ...)
      ip2_f <- ip2_est$ip_f
      ip2 <- ip2_est$ip
      
    } else {
      
      stop("Invalid 'input.type' value!")
      
    } 
    
    #############
    ## Mapping ##
    #############
    
    ## Homography Transformation Method
    cat("\nMapping d2 ideal points on d1 ideal points space...\n\n")
    if (bridge.method=="procrustes") {
      
      bridge.model <- bridge.procrustes(ip1=ip1,
                                        ip2=ip2,
                                        anchorrows.ip1=anchorrows.ip1,
                                        anchorrows.ip2=anchorrows.ip2, 
                                        opt = tr.opt,
                                        opt.iter.n = tr.opt.iter.n,
                                        opt.sample.n = tr.opt.sample.n,
                                        opt.th.inline = tr.opt.th.inline,
                                        blend = tr.blend,
                                        blend.th1 = tr.blend.th1, 
                                        blend.th2 = tr.blend.th2)
      ip2_trans <- bridge.model$ip2_trans
      
    }
    if (bridge.method=="homography") {
      
      bridge.model <- bridge.homography(ip1=ip1,
                                        ip2=ip2,
                                        anchorrows.ip1=anchorrows.ip1,
                                        anchorrows.ip2=anchorrows.ip2, 
                                        opt = tr.opt,
                                        opt.iter.n = tr.opt.iter.n,
                                        opt.sample.n = tr.opt.sample.n,
                                        opt.th.inline = tr.opt.th.inline,
                                        blend = tr.blend,
                                        blend.th1 = tr.blend.th1, 
                                        blend.th2 = tr.blend.th2)
      ip2_trans <- bridge.model$ip2_trans
      
    }
    if (bridge.method=="olsmap") {
      
      bridge.model <- bridge.olsmap(ip1=ip1,
                                        ip2=ip2,
                                        anchorrows.ip1=anchorrows.ip1,
                                        anchorrows.ip2=anchorrows.ip2, 
                                        opt = tr.opt,
                                        opt.iter.n = tr.opt.iter.n,
                                        opt.sample.n = tr.opt.sample.n,
                                        opt.th.inline = tr.opt.th.inline,
                                        blend = tr.blend,
                                        blend.th1 = tr.blend.th1, 
                                        blend.th2 = tr.blend.th2)
      ip2_trans <- bridge.model$ip2_trans
      
    }
    
    # Update respdt for anchor rows
    respdt$isanchor[anchorrows.ip1[which(anchorrows.ip1 <= nrow(d1))]] <- 1
    respdt$isanchor[anchorrows.ip2[which(anchorrows.ip2 <= nrow(d2))]+nrow(d1)] <- 1

    ## Add IP values to respdt
    for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")] <- NA
    for(i in 1:ncol(ip1)) respdt[,paste0("ip1_coord",i,"D")] <- NA
    for(i in 1:ncol(ip1)) respdt[,paste0("ip2_coord",i,"D")] <- NA
    if (anchors.method=="selectrows") {
      
      for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")][which(respdt$data==2)] <- ip2_trans[,i]
      for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")][which(respdt$data==1)] <- ip1[,i]
      for(i in 1:ncol(ip1)) respdt[,paste0("ip1_coord",i,"D")][which(respdt$data==1)] <- ip1[,i]
      for(i in 1:ncol(ip1)) respdt[,paste0("ip2_coord",i,"D")][which(respdt$data==2)] <- ip2[,i]
      # If anchors.method=="selectrows", anchors are duplicated in d2. 
      # Therefore bridged ideal points of anchors in d2 are replaced with NA.
      for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")][which(respdt$data==2 & respdt$isanchor==1)] <- NA
      
    } else {
      
      for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")][which(respdt$data==2|respdt$isanchor==1)] <- ip2_trans[,i]
      for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")][which(respdt$data==1|respdt$isanchor==1)] <- ip1[,i]
      for(i in 1:ncol(ip1)) respdt[,paste0("ip1_coord",i,"D")][which(respdt$data==1|respdt$isanchor==1)] <- ip1[,i]
      for(i in 1:ncol(ip1)) respdt[,paste0("ip2_coord",i,"D")][which(respdt$data==2|respdt$isanchor==1)] <- ip2[,i]
      
    }

    cat("DONE!\n\n")
    ## Compile Output
    out <- list(bridge.data = respdt, 
                bridge.model = bridge.model,
                bridge.method = bridge.method, 
                ip.model.d1 = ip1_f,
                ip.model.d2 = ip2_f,
                ip.method = ip.method,
                ip.dims = ip.dims,
                input.type = input.type,
                anchors.method = anchors.method)
    return(out)
  
  # Mapping using Model Prediction from Only One Group (Only Bayesian IRT)        
  } else if (bridge.method=="modelmap") {
    
    if (input.type=="idealpoints"|ip.method!="irtMCMC") {
      stop("To use bridge.method='modelmap', input.type must be 'responses' and ip.method must be 'irtMCMC'!")
    }
    
    ip_mmap_f <- ideal(rollcall(rbind(d1,d2)), d=ip.dims, 
                       use.voter=c(rep(TRUE,nrow(d1)),rep(FALSE,nrow(d2))), ...)
    ip_mmap <- ip_mmap_f$xbar
    
    ## Add IP values to respdt
    for(i in 1:ncol(ip_mmap)) respdt[,paste0("bridged",i,"D")] <- ip_mmap[,i]
    for(i in 1:ncol(ip_mmap)) respdt[,paste0("ip1_coord",i,"D")] <- NA
    for(i in 1:ncol(ip_mmap)) respdt[,paste0("ip2_coord",i,"D")] <- NA
    
    cat("DONE!\n\n")
    
    ## Compile Output
    out <- list(bridge.data = respdt, 
                bridge.model = ip_mmap_f,
                bridge.method = bridge.method, 
                ip.model.d1 = NULL,
                ip.model.d2 = NULL,
                ip.method = ip.method,
                ip.dims = ip.dims,
                input.type = input.type,
                anchors.method = NULL)
    return(out)

  } else if (bridge.method=="anchoredpooling") {
    
    if (ncol(d1)!=ncol(d2)) stop("d1 and d2 must have the same number of columns!")
    
    if (input.type=="idealpoints") {
      
      stop('input.type=="idealpoints" is incompatible with bridge.method=="joint"!')
      
    } else if (input.type=="responses") {
      
      if (anchors.method!="selectrows") stop("anchors.method must be 'selectrows'!")
      
      if (length(anchors.selectrows.d1)!=length(anchors.selectrows.d2)) {
        stop("anchors.selectrows.d1 and anchors.selectrows.d2 must have the same length!")
      }
      
      # Update respdt
      respdt$isanchor[anchors.selectrows.d1] <- 1
      respdt$isanchor[anchors.selectrows.d2+nrow(d1)] <- 1
      
      ## Dataset for Ideal Point Computation
      d1x <- cbind(d1[-anchors.selectrows.d1,], 
                   matrix(NA,ncol=ncol(d2),nrow=nrow(d1[-anchors.selectrows.d1,])))
      d2x <- cbind(matrix(NA,ncol=ncol(d1),nrow=nrow(d2[-anchors.selectrows.d2,])),
                   d2[-anchors.selectrows.d2,])
      dx <- rbind(d1x,d2x,cbind(d1[anchors.selectrows.d1,],
                                d2[anchors.selectrows.d2,]))

      cat("Generating bridged ideal points...\n")
      
      # Anchor Identifier (NA in joint)
      respdt$isanchor <- NA
      
      # Estimate Ideal Points
      ip_pooled_est <- ipest(dx, method=ip.method, 
                            dims=ip.dims, polarity=ip.polarity.d1, ...)
      ip_pooled_f <- ip_pooled_est$ip_f
      ip_pooled <- ip_pooled_est$ip
      
    } else {
      
      stop("Invalid 'input.type' value!")
      
    }
    
    ## Add IP values to respdt
    acrows <- c(anchors.selectrows.d1,anchors.selectrows.d2+nrow(d1))
    for(i in 1:ncol(ip_pooled)) {
      respdt[,paste0("bridged",i,"D")] <- NA
      respdt[-acrows,paste0("bridged",i,"D")] <- 
        ip_pooled[seq(1,nrow(d1x)+nrow(d2x)),i]
      respdt[anchors.selectrows.d1,paste0("bridged",i,"D")] <- 
        ip_pooled[-seq(1,nrow(d1x)+nrow(d2x)),i]
    }
    for(i in 1:ncol(ip_pooled)) respdt[,paste0("ip1_coord",i,"D")] <- NA
    for(i in 1:ncol(ip_pooled)) respdt[,paste0("ip2_coord",i,"D")] <- NA
    
    cat("DONE!\n\n")
    
    ## Compile Output
    out <- list(bridge.data = respdt, 
                bridge.model = ip_pooled_f,
                bridge.method = bridge.method, 
                ip.model.d1 = NULL,
                ip.model.d2 = NULL,
                ip.method = ip.method,
                ip.dims = ip.dims,
                input.type = input.type,
                anchors.method = "selectrows")
    return(out)
    
  } else {
    
    stop("Invalid 'bridge.method' value!")
    
  }
  
  
  
}

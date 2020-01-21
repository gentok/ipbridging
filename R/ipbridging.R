#' Bridging Ideal Point Estimates from Two Data Sets
#' 
#' @description Toolbox function that Bridges ideal point estimates computed from 
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
#'   \item \code{"homography"} (default): Homography transformation method. 
#'   Non-parametric method to find optimal transformation matrix to bridge ideal point 
#'   estimates (calls \code{\link{bridge.homography}}).
#'   \item \code{"olsmap"}: Use OLs regression to map \code{d2} ideal point coordinates 
#'   on \code{d1} ideal point space. The procedure described in Shor et al. (2010). 
#'   \item \code{"pooling"}: Pooling method. Simple method to combine \code{d1} and 
#'   \code{d2} to jointly compute ideal points. The procedure described in Jessee (2016).
#' }
#' @param anchors.method Method to select anchors to bridge ideal points. 
#' Following values are currently available.
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
#'   \item \code{"random"} (default): randomly sample subset of dataset and 
#'   use them as anchors. The sample size is determined by \code{anchors.subsample.pr}.
#'   \item \code{"selectrows"}: Use this option if specific anchors are already known (
#'   but not shared between datasets). Row numbers for anchors must be set 
#'   in \code{anchors.selectrows.d1} and \code{anchors.selectrows.d2}, respectively
#'   (Two vectors does not have to have the same length).
#' }
#' @param anchors.subsample.pr Proportion of dataset (i.e., \code{d1} and \code{d2}) 
#' to be sampled as anchor. The propotion is determined by the dataset of smaller size.
#' Used when \code{anchors.method=="subsample"} and \code{anchors.subsample.method=="random"}.
#' @param anchors.selectrows.d1 Must be positive integers.  
#' Specify the row numbers of anchors in \code{d1}. 
#' If \code{anchors.method=="selectrows"}, must be the same length and with perfect correspondence with \code{anchors.selectrows.d2}.
#' If \code{anchors.method=="subsample"} and \code{anchors.subsample.method=="selectrows"}, 
#' the above condition is not required. 
#' @param anchors.selectrows.d2 Must be positive integers.  
#' Specify the row numbers of anchors in \code{d2}.
#' @param anchors.newdata Matrix or data.frame of anchors. Must have identical 
#' columns as \code{d1} and \code{d2}. Used if \code{anchors.method=="newdata"}. 
#' @param input.type The input type of \code{d1} and \code{d2}. If \code{"responses"} (default), 
#' dataset values are responses to multiple issues or votes. 
#' If \code{"idealpoints"}, dataset values are ideal points.   
#' @param ip.method Method of ideal point computation (used only when \code{input.type=="response"}. Following values are currently available.
#' \itemize{
#'   \item \code{"ooc"} (default): Ordered optimal classification. Calls \code{\link{oocflex}} function.
#'   \item \code{"oc"}: Optimal classification. Calls \code{\link[oc]{oc}} function from \code{oc} package.
#'   \item \code{"wnominate"}: W-NOMINATE. Calls \code{\link[wnominate]{wnominate}} function from \code{wnominate} package.
#'   \item \code{"irtMCMC"}: IRT model via Markov chain Monte Carlo method. Calls \code{\link[pscl]{ideal}} function from \code{pscl} package.
#' }
#' For \code{ip.method=="ooc"}, the inputs \code{d1} and \code{d2} must be positive integer. 
#' For all other methods, the inputs must be a dummy variable indication 0=Nay and 1=Yea (and 9 for unit nonresponse).
#' All methods allow item non-response as NA.
#' @param ip.dims Number of dimension in ideal point computation. Must be a positive 
#' integer between 1 and 10. If \code{bridge.method=="homography"}, only 2 is accepted 
#' at this point.  
#' @param hg.opt.iter.n Used if \code{bridge.method=="homography"}. 
#' Number of iteration in the optimization of transformation matrix.
#' @param hg.opt.sample.n Used if \code{bridge.method=="homography"}. 
#' Size of anchoring respondents to be sub-sampled at each iteration of 
#' optimization.
#' @param hg.opt.th.inline Used if \code{bridge.method=="homography"}. 
#' Upper bound to determine inline respondents at each iteration of 
#' optimization. A respondent is considered 'inline' if the distance between transformed
#' \code{ip2} and \code{ip1} goes below this threshold.
#' @param hg.blend.th1 Used if \code{bridge.method=="homography"}. 
#' First threshold used in 'blending' procedure. If minimum difference 
#' between transformed \code{ip2} and \code{ip1} goes below this threshold, the 
#' transformed \code{ip2} is replaced with the closest \code{ip1}. 
#' @param hg.blend.th2 Used if \code{bridge.method=="homography"}. Second threshold used in 'blending' procedure. If minimum difference 
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
#' @seealso \code{\link{bridge.homography}}, \code{\link{oocflex}}, \code{\link[oc]{oc}}, 
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
#' 
#' @export

ipbridging <- function(d1, d2, 
                       bridge.method = "homography", 
                       anchors.method = "subsample",
                       anchors.subsample.pr = 0.1,
                       anchors.subsample.method = "random",
                       anchors.selectrows.d1 = 1:100,
                       anchors.selectrows.d2 = 1:100,
                       anchors.newdata = NULL,
                       input.type="responses",
                       ip.method="ooc", 
                       ip.dims=2,
                       hg.opt.iter.n = 10000,
                       hg.opt.sample.n = 30,
                       hg.opt.th.inline = 0.5,
                       hg.blend.th1 = 0.05, 
                       hg.blend.th2 = 0.15,
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
      
  ## Pooling As a Bridging Method
  if (bridge.method=="pooling") {

    if (ncol(d1)!=ncol(d2)) stop("d1 and d2 must have the same number of columns!")
    
    if (input.type=="idealpoints") {
      
      stop('input.type=="idealpoints" is incompatible with bridge.method=="pooling"!')
      
    } else if (input.type=="responses") {
      
      cat("Generating bridged ideal points...\n")
      
      # Anchor Identifier (NA in pooling)
      respdt$isanchor <- NA

      # Estimate Ideal Points
      if (ip.method=="ooc") {
        
        ip_pooled_f <- oocflex(rbind(d1,d2), dims=ip.dims, ...)
        ip_pooled <- ip_pooled_f$respondents[,grepl("coord", colnames(ip_pooled_f$respondents))]
        
      } else if (ip.method=="oc") {
        
        ip_pooled_f <- oc(rollcall(rbind(d1,d2)), dims=ip.dims, ...)
        ip_pooled <- ip_pooled_f$legislators[,grepl("coord", colnames(ip_pooled_f$legislators))]

      } else if (ip.method=="wnominate") {
        
        ip_pooled_f <- wnominate(rollcall(rbind(d1,d2)), dims=ip.dims, ...)
        ip_pooled <- ip_pooled_f$legislators[,grepl("coord", colnames(ip_pooled_f$legislators))]
        
      } else if (ip.method=="irtMCMC") {
        
        ip_pooled_f <- ideal(rollcall(rbind(d1,d2)), d=ip.dims, ...)
        ip_pooled <- ip_pooled_f$xbar
        
      } else {
        
        stop("Invalid 'ip.method' value!")
        
      }

    } else {
      
      stop("Invalid 'input.type' value!")
      
    }
    
    ## Add IP values to respdt
    for(i in 1:ncol(ip_pooled)) respdt[,paste0("bridged",i,"D")] <- ip_pooled[,i]
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
                anchors.method = NULL)
    return(out)
    
  ## All Other Bridging Methods using Linear Mapping with Anchors    
  } else if (bridge.method%in%c("homography","olsmap")) {
    
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

        # Update respdt
        respdt$isanchor[anchors.selectrows.d1] <- 1
        respdt$isanchor[anchors.selectrows.d2+nrow(d1)] <- 1
        
      }
      
      
    } else if (input.type=="responses") {
      
      ## Set Anchors
      cat("Setting anchors...\n\n")
      if (anchors.method=="newdata") {
        
        if (ncol(d1)!=ncol(d2)) stop("d1 and d2 must have the same number of columns!")
        
        if (is.null(anchors.newdata)) {
          
          stop("anchors.newdata is NULL!")
          
        } else if (!class(anchors.newdata)[1]%in%c("matrix","data.frame")) {
          
          stop("anchors.newdata must be matrix or data.frame!")
          
        } else if (ncol(anchors.newdata)!=ncol(d1)) {
          
          stop("anchors.newdata must have same number of column as ")
          
        }
        if (is.data.frame(anchors.newdata)) anchors.newdata <- as.matrix(anchors.newdata) 
        
        # Set Anchors data
        ac <- anchors.newdata
        
        # Additional rows in respdt
        respdt_add <- data.frame(allid = seq(nrow(d1)+nrow(d2)+1,nrow(d1)+nrow(d2)+nrow(ac)),
                                 subid = seq(1,nrow(ac)),
                                 data = 3,
                                 isanchor = 1)
        respdt <- rbind(respdt, respdt_add)
        
        ## Dataset for Ideal Point Computation
        d1x <- rbind(d1, ac)
        d2x <- rbind(d2, ac)
        anchorrows.ip1 <- seq(nrow(d1)+1,nrow(d1x))
        anchorrows.ip2 <- seq(nrow(d2)+1,nrow(d2x))
        
      } else if (anchors.method=="selectrows") { 
        
        if (length(anchors.selectrows.d1)!=length(anchors.selectrows.d2)) {
          stop("anchors.selectrows.d1 and anchors.selectrows.d2 must have the same length!")
        }
        
        # Update respdt
        respdt$isanchor[anchors.selectrows.d1] <- 1
        respdt$isanchor[anchors.selectrows.d2+nrow(d1)] <- 1
        
        ## Dataset for Ideal Point Computation
        d1x <- d1
        d2x <- d2
        anchorrows.ip1 <- anchors.selectrows.d1
        anchorrows.ip2 <- anchors.selectrows.d2
        
      } else if (anchors.method=="subsample") {
        
        ## Sample subset of rows and use it as anchor
        if (anchors.subsample.method=="random") {
          
          ## Random Sample of Rows
          sample.size <-  floor(min(nrow(d1),nrow(d2))*anchors.subsample.pr)
          anchors.samplerows.d1 <- sample.int(nrow(d1), size = sample.size)
          anchors.samplerows.d2 <- sample.int(nrow(d2), size = sample.size)
          
        } else if (anchors.subsample.method=="selectrows") { 
          
          ## Subsample by manually selecting rows
          anchors.samplerows.d1 <- anchors.selectrows.d1
          anchors.samplerows.d2 <- anchors.selectrows.d2

        } else {
          
          stop("Invalid 'anchors.subsample.method' value!")
          
        }
        
        # Update respdt
        respdt$isanchor[anchors.samplerows.d1] <- 1
        respdt$isanchor[anchors.samplerows.d2+nrow(d1)] <- 1
        
        ## Update Dataset for Ideal Point Computation
        d1x <- rbind(d1, d2[anchors.samplerows.d2,])
        d2x <- rbind(d2, d1[anchors.samplerows.d1,])
        anchorrows.ip1 <- c(anchors.samplerows.d1, seq(nrow(d1)+1,nrow(d1x)))
        anchorrows.ip2 <- c(anchors.samplerows.d2, seq(nrow(d2)+1,nrow(d2x)))
        
      } else {
        
        stop("Invalid 'anchors.method' value!")
        
      }
      
      ## Compute Ideal Points
      cat("Generating ideal points on subsets...\n")
      if (ip.method=="ooc") {
        
        ip1_f <- oocflex(d1x, dims=ip.dims, ...)
        ip2_f <- oocflex(d2x, dims=ip.dims, ...)
        ip1 <- ip1_f$respondents[,grepl("coord", colnames(ip1_f$respondents))]
        ip2 <- ip2_f$respondents[,grepl("coord", colnames(ip2_f$respondents))]
        
      } else if (ip.method=="oc") {
        
        ip1_f <- oc(rollcall(d1x), dims=ip.dims, ...)
        ip2_f <- oc(rollcall(d2x), dims=ip.dims, ...)
        ip1 <- ip1_f$legislators[,grepl("coord", colnames(ip1_f$legislators))]
        ip2 <- ip2_f$legislators[,grepl("coord", colnames(ip2_f$legislators))]
        
      } else if (ip.method=="wnominate") {
        
        ip1_f <- wnominate(rollcall(d1x), dims=ip.dims, ...)
        ip2_f <- wnominate(rollcall(d2x), dims=ip.dims, ...)
        ip1 <- ip1_f$legislators[,grepl("coord", colnames(ip1_f$legislators))]
        ip2 <- ip2_f$legislators[,grepl("coord", colnames(ip2_f$legislators))]
        
      } else if (ip.method=="irtMCMC") {
        
        ip1_f <- ideal(rollcall(d1x), d=ip.dims, ...)
        ip2_f <- ideal(rollcall(d2x), d=ip.dims, ...)
        ip1 <- ip1_f$xbar
        ip2 <- ip2_f$xbar
        
      } else {
        
        stop("Invalid 'ip.method' value!")
        
      }
      
      #############
      ## Mapping ##
      #############
      
      ## Homography Transformation Method
      cat("\nMapping d2 ideal points on d1 ideal points space...\n\n")
      if (bridge.method=="homography") {
        
        bridge.model <- bridge.homography(ip1=ip1,
                                    ip2=ip2,
                                    anchorrows.ip1=anchorrows.ip1,
                                    anchorrows.ip2=anchorrows.ip2, 
                                    opt.iter.n = hg.opt.iter.n,
                                    opt.sample.n = hg.opt.sample.n,
                                    opt.th.inline = hg.opt.th.inline,
                                    blend.th1 = hg.blend.th1, 
                                    blend.th2 = hg.blend.th2)
        ip2_trans <- bridge.model$ip2_trans
        
      }
      if (bridge.method=="olsmap") {
        
        bridge.model <- bridge.olsmap(ip1=ip1,
                                      ip2=ip2,
                                      anchorrows.ip1=anchorrows.ip1,
                                      anchorrows.ip2=anchorrows.ip2)
        ip2_trans <- bridge.model$ip2_trans
        
      }

      ## Add IP values to respdt
      for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")] <- NA
      for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")][which(respdt$data==2|respdt$isanchor==1)] <- ip2_trans[,i]
      for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")][which(respdt$data==1|respdt$isanchor==1)] <- ip1[,i]
      if (anchors.method=="selectrows") {
        # If anchors.method=="selectrows", anchors are duplicated in d2. 
        # Therefore bridged ideal points of anchors in d2 are replaced with NA.
        for(i in 1:ncol(ip1)) respdt[,paste0("bridged",i,"D")][which(respdt$data==2 & respdt$isanchor==1)] <- NA
      }
      for(i in 1:ncol(ip1)) respdt[,paste0("ip1_coord",i,"D")] <- NA
      for(i in 1:ncol(ip1)) respdt[,paste0("ip1_coord",i,"D")][which(respdt$data==1|respdt$isanchor==1)] <- ip1[,i]
      for(i in 1:ncol(ip1)) respdt[,paste0("ip2_coord",i,"D")] <- NA
      for(i in 1:ncol(ip1)) respdt[,paste0("ip2_coord",i,"D")][which(respdt$data==2|respdt$isanchor==1)] <- ip2[,i]
      
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

    } else {
      
      stop("Invalid 'input.type' value!")
      
    }    
        
  } else {
    
    stop("Invalid 'bridge.method' value!")
    
  }
  
  
  
}

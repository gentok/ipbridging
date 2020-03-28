#' Setting Anchors to Combine Two Datasets
#' 
#' @description Toolbox function that sets anchors to bridge two datasets.
#' 
#' @param d1 Matrix or data.frame of 'reference' responses.
#' (i.e., In mapping methods, ideal points based on \code{d2} will be 
#' transformed and mapped to ideal point space based on \code{d1}).
#' Rows are resepondents and columns are questions/issues. 
#' @param d2 Matrix or data.frame of responses to be transformed.
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
#' @param anchors.subsample.pr Numeric vector of length 1 or 2. 
#' Used when  \code{anchors.method=="subsample"}. 
#' Proportion of dataset (i.e., \code{d1} and \code{d2}) to be sampled as anchor. 
#' When length is 1, the propotion is determined by the dataset of smaller size if 
#' \code{anchors.subsample.method=="random"}, specific dataset's size if 
#' \code{anchors.subsample.method} is \code{"random.d1"} or \code{"random.d2"}. 
#' When length is 2 and \code{anchors.subsample.method=="random"}, the first 
#' element determines the proportion in \code{d1}, the second element determines 
#' the proportion in \code{d2}.
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
#' @param anchors.selectrows.data Required if \code{anchors.method=="selectrows"}. 
#' Specify the dataset where anchors are originally found. 
#' Must be the same length as \code{anchors.selectrows.d1}. Allowed values are 1 (=\code{d1}), 
#' 2 (=\code{d2}) and 3 (=new data).  
#' 
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{bridge.data}: data.frame to store the ideal point estimates in the next method.
#'   \item \code{d.ip1}: \code{d1} with anchors.
#'   \item \code{d.ip2}: \code{d2} with anchors.
#'   \item \code{anchorrows.ip1}: Anchor locations in \code{d.ip1}.
#'   \item \code{anchorrows.ip2}: Anchor locations in \code{d.ip2}.
#'   \item \code{anchorrows.ip1.data}: Original dataset for anchors identified in \code{anchorrows.ip1}.
#'   \item \code{anchorrows.ip2.data}: Original dataset for anchors identified in \code{anchorrows.ip2}.
#' }
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' 
#' @seealso \code{\link{ipbridging}}, \code{\link{ipest}}
#' 
#' @import stats
#' 
#' @export

setanchors <- function(d1,d2,
                       anchors.method = "subsample",
                       anchors.subsample.method = "random",
                       anchors.subsample.pr = 0.1,
                       anchors.subsample.wgt.d1 = NULL,
                       anchors.subsample.wgt.d2 = NULL,
                       anchors.selectrows.d1 = 1:100,
                       anchors.selectrows.d2 = 1:100,
                       anchors.selectrows.data = NULL,
                       anchors.newdata = NULL
                       ) {
  
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
  
  ## Set Anchors
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
    
    # Anchor ID
    respdt$anchorid <- NA
    respdt$anchorid[which(repdt$isanchor==1)] <- 
      seq(1,length(which(respdt$isanchor==1))) 
    
    ## Dataset for Ideal Point Computation
    d1x <- rbind(d1, ac)
    d2x <- rbind(d2, ac)
    
    # Row Number for Anchors
    anchorrows.ip1 <- seq(nrow(d1)+1,nrow(d1x))
    anchorrows.ip2 <- seq(nrow(d2)+1,nrow(d2x))
    
    ## Original Dataset of Anchors 
    anchorrows.ip1.data <- rep(3, length(anchorrows.ip1))
    anchorrows.ip2.data <- rep(3, length(anchorrows.ip2))
    
  } else if (anchors.method=="selectrows") { 
    
    if (length(anchors.selectrows.d1)!=length(anchors.selectrows.d2)) {
      stop("anchors.selectrows.d1 and anchors.selectrows.d2 must have the same length!")
    }
    
    if (is.null(anchors.selectrows.data)) {
      stop("anchors.selectrows.data must be set!")
    }
    
    if (length(anchors.selectrows.data)!=length(anchors.selectrows.d1)) {
      stop("anchors.selectrows.data and anchors.selectrows.d1 must have the same length!")
    }
    
    # Respondents Data
    # d1
    respdt1 <- data.frame(allid = NA, subid = NA,
                          data = rep(1,nrow(d1)), 
                          isanchor = 0, anchorid = NA)
    respdt1$isanchor[anchors.selectrows.d1] <- 1
    for(i in 1:length(anchors.selectrows.d1)) {
      respdt1$anchorid[anchors.selectrows.d1[i]] <- i
      respdt1$data[anchors.selectrows.d1[i]] <- anchors.selectrows.data[i]
    }
    respdt1 <- respdt1[which(respdt1$data==1),]
    # d2
    respdt2 <- data.frame(allid = NA, subid = NA,
                          data = rep(2,nrow(d2)), isanchor = 0)
    respdt2$isanchor[anchors.selectrows.d2] <- 1
    for(i in 1:length(anchors.selectrows.d2)) {
      respdt2$anchorid[anchors.selectrows.d2[i]] <- i
      respdt2$data[anchors.selectrows.d2[i]] <- anchors.selectrows.data[i]
    }
    respdt2 <- respdt2[which(respdt2$data==2),]
    ## Combine d1 and d2
    respdt <- rbind(respdt1, respdt2)
    # d3 (if any)
    if(length(which(anchors.selectrows.data==3))>0) {
      # Additional rows in respdt
      respdt_add <- data.frame(allid = NA,
                               subid = seq(1,length(which(anchors.selectrows.data==3))),
                               data = 3,
                               isanchor = 1,
                               anchorid = which(anchors.selectrows.data==3))
      respdt <- rbind(respdt, respdt_add)
    }
    
    ## Dataset for Ideal Point Computation
    d1x <- rbind(d1[which(!seq(1,nrow(d1)) %in% anchors.selectrows.d1[anchors.selectrows.data!=1])],
                 d1[anchors.selectrows.d1[which(anchors.selectrows.data==2)]],
                 d1[anchors.selectrows.d1[which(anchors.selectrows.data==3)]])
    d2x <- rbind(d2[which(!seq(1,nrow(d2)) %in% anchors.selectrows.d2[anchors.selectrows.data!=2])],
                 d2[anchors.selectrows.d2[which(anchors.selectrows.data==1)]],
                 d2[anchors.selectrows.d2[which(anchors.selectrows.data==3)]])
    
    ## Anchor Rows
    ### Data 1
    anchorrows.ip1 <- as.numeric(sapply(which(anchors.selectrows.data==1), function(k) which(respdt[respdt$data==1,]$anchorid==k)))
    if (length(which(anchors.selectrows.data==2))>0) {
      anchorrows.ip1 <- c(anchorrows.ip1,
                          seq(nrow(d1)-length(which(anchors.selectrows.data!=1))+1,
                              nrow(d1)-length(which(anchors.selectrows.data!=1))+length(which(anchors.selectrows.data==2))))
    }
    if (length(which(anchors.selectrows.data==3))>0) {
      anchorrows.ip1 <- c(anchorrows.ip1,
                          seq(nrow(d1)-length(which(anchors.selectrows.data==3))+1,
                              nrow(d1)-length(which(anchors.selectrows.data==3))+length(which(anchors.selectrows.data==3))))
    }
    ### Data 2
    anchorrows.ip2 <- as.numeric(sapply(which(anchors.selectrows.data==2), function(k) which(respdt[respdt$data==2,]$anchorid==k)))
    if (length(which(anchors.selectrows.data==1))>0) {
      anchorrows.ip2 <- c(anchorrows.ip2,
                          seq(nrow(d1)-length(which(anchors.selectrows.data!=2))+1,
                              nrow(d1)-length(which(anchors.selectrows.data!=2))+length(which(anchors.selectrows.data==1))))
    }
    if (length(which(anchors.selectrows.data==3))>0) {
      anchorrows.ip2 <- c(anchorrows.ip2,
                          seq(nrow(d1)-length(which(anchors.selectrows.data==3))+1,
                              nrow(d1)-length(which(anchors.selectrows.data==3))+length(which(anchors.selectrows.data==3))))
    }
    
    ## Original Dataset of Anchors 
    anchorrows.ip1.data <- c(rep(1, length(which(anchors.selectrows.data==1))),
                             rep(2, length(which(anchors.selectrows.data==2))),
                             rep(3, length(which(anchors.selectrows.data==3))))
    anchorrows.ip2.data <- c(rep(2, length(which(anchors.selectrows.data==2))),
                             rep(1, length(which(anchors.selectrows.data==1))),
                             rep(3, length(which(anchors.selectrows.data==3))))
    
  } else if (anchors.method=="subsample") {
    
    # Respondents Data
    respdt <- data.frame(allid = seq(1,nrow(d1)+nrow(d2)),
                         subid = c(seq(1,nrow(d1)),seq(1,nrow(d2))),
                         data = c(rep(1, nrow(d1)),rep(2, nrow(d2))),
                         isanchor = 0)
    
    ## Sample subset of rows and use it as anchor
    if (anchors.subsample.method=="random") {
      
      ## Random Sample of Rows
      if (length(anchors.subsample.pr)>2) stop("Length of anchors.subsample.pr must be <=2.")
      if (length(anchors.subsample.pr)==1) {
        anchors.sample.size <-  rep(floor(min(nrow(d1),nrow(d2))*anchors.subsample.pr),2)
      } else {
        anchors.sample.size <-  c(nrow(d1)*anchors.subsample.pr[1],
                                  nrow(d2)*anchors.subsample.pr[2])
      }
      anchors.samplerows.d1 <- sample.int(nrow(d1), 
                                          size = anchors.sample.size[1],
                                          prob = anchors.subsample.wgt.d1)
      anchors.samplerows.d2 <- sample.int(nrow(d2), 
                                          size = anchors.sample.size[2],
                                          prob = anchors.subsample.wgt.d2)
      
    } else if (anchors.subsample.method=="random.d1") {
      
      ## Random Sample of Rows
      anchors.sample.size <-  floor(nrow(d1)*anchors.subsample.pr)
      anchors.samplerows.d1 <- sample.int(nrow(d1), 
                                          size = anchors.sample.size,
                                          prob = anchors.subsample.wgt.d1)
      anchors.samplerows.d2 <- integer() 
      
    } else if (anchors.subsample.method=="random.d2") {
      
      ## Random Sample of Rows
      anchors.sample.size <-  floor(nrow(d2)*anchors.subsample.pr)
      anchors.samplerows.d1 <- integer() 
      anchors.samplerows.d2 <- sample.int(nrow(d2), 
                                          size = anchors.sample.size,
                                          prob = anchors.subsample.wgt.d2)
      
    } else if (anchors.subsample.method=="selectrows") { 
      
      ## Subsample by manually selecting rows
      anchors.samplerows.d1 <- anchors.selectrows.d1
      anchors.samplerows.d2 <- anchors.selectrows.d2
      
    } else {
      
      stop("Invalid 'anchors.subsample.method' value!")
      
    }
    
    # Update respdt
    respdt$isanchor[respdt$data==1][anchors.samplerows.d1] <- 1
    respdt$isanchor[respdt$data==2][anchors.samplerows.d2] <- 1
    
    ## Update Dataset for Ideal Point Computation
    d1x <- rbind(d1, d2[anchors.samplerows.d2,])
    d2x <- rbind(d2, d1[anchors.samplerows.d1,])
    if (nrow(d1x)>nrow(d1)) {
      anchorrows.ip1 <- c(anchors.samplerows.d1, seq(nrow(d1)+1,nrow(d1x)))
      anchorrows.ip1.data <- c(rep(1, length(anchors.samplerows.d1)),
                               rep(2, length(seq(nrow(d1)+1,nrow(d1x)))))
    } else {
      anchorrows.ip1 <- anchors.samplerows.d1
      anchorrows.ip1.data <- c(rep(1, length(anchors.samplerows.d1)))
    }
    if (nrow(d2x)>nrow(d2)) {
      anchorrows.ip2 <- c(anchors.samplerows.d2, seq(nrow(d2)+1,nrow(d2x)))
      anchorrows.ip2.data <- c(rep(2, length(anchors.samplerows.d2)),
                               rep(1, length(seq(nrow(d2)+1,nrow(d2x)))))
    } else {
      anchorrows.ip2 <- anchors.samplerows.d2
      anchorrows.ip2.data <- c(rep(2, length(anchors.samplerows.d2)))
    }
    
  } else {
    
    stop("Invalid 'anchors.method' value!")
    
  }
  
  out <- list(bridge.data=respdt, 
              d.ip1=d1x, d.ip2=d2x,
              anchorrows.ip1=anchorrows.ip1,
              anchorrows.ip2=anchorrows.ip2,
              anchorrows.ip1.data=anchorrows.ip1.data,
              anchorrows.ip2.data=anchorrows.ip2.data,
              anchors.method = anchors.method)
  
  if (anchors.method=="subsample") {
    out$anchors.subsample.method = anchors.subsample.method
    out$anchors.subsample.pr = anchors.subsample.pr
    out$anchors.subsample.wgt.d1 = anchors.subsample.wgt.d1
    out$anchors.subsample.wgt.d2 = anchors.subsample.wgt.d1
    if (anchors.subsample.method=="selectrows") {
      out$anchors.selectrows.d1 = anchors.selectrows.d1
      out$anchors.selectrows.d2 = anchors.selectrows.d2
    }
  } else if (anchors.method=="selectrows") {
    out$anchors.selectrows.d1 = anchors.selectrows.d1
    out$anchors.selectrows.d2 = anchors.selectrows.d2
    out$anchors.selectrows.data = anchors.selectrows.data
  } else if (anchors.method=="newdata") {
    out$anchors.newdata = anchors.newdata
  }
  
  class(out) <- c("anchors.out", class(out))
  return(out)
  
}
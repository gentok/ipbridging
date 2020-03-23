#' Estimating Ideal Points
#' 
#' @description Toolbox function that estimates ideal points using different methods.
#' 
#' @param d A matrix with numeric values. With following conditions:
#' \itemize{
#'   \item Missing values (if any) must be set to \code{NA}.
#'   \item If variables are binary, values should be conisted from 0 and 1.
#'   \item If variables are ordered, values should be positive integers starting from 1.
#' }
#' 
#' @param method Method of ideal point computation (used only when \code{input.type=="response"}. Following values are currently available.
#' \itemize{
#'   \item \code{"ooc"} (default): Ordered optimal classification. Calls \code{\link{oocflex}} function.
#'   \item \code{"blackbox"}: Blackbox scaling. Calls \code{\link[basicspace]{blackbox}} function from \code{blackbox} package (\code{minscale=10} by default).
#'   \item \code{"oc"}: Optimal classification. Calls \code{\link[oc]{oc}} function from \code{oc} package (\code{polarity=c(1,1)} by default). 
#'   \item \code{"wnominate"}: W-NOMINATE. Calls \code{\link[wnominate]{wnominate}} function from \code{wnominate} package  (\code{polarity=c(1,1)} by default).
#'   \item \code{"irtMCMC"}: IRT model via Markov chain Monte Carlo method. 
#'   Calls \code{\link[pscl]{ideal}} function from \code{pscl} package if \code{irt.type=="binary"};
#'   Calls \code{\link[MCMCpack]{MCMCfactanal}} function from \code{MCMCpack} package if \code{irt.type=="ordered"} (still provisional).
#'   \item \code{"irtEM"}: IRT model via EM algorithm method (only one dimensional ideal points. This option is still provisional).
#'   Calls \code{\link[emIRT]{binIRT}} function from \code{emIRT} package if \code{irt.type=="binary"};
#'   Calls \code{\link[emIRT]{ordIRT}} function from \code{emIRT} package if \code{irt.type=="ordered"} (Only 3 category variable can be used).
#' }
#' For \code{method=="ooc"}, the inputs \code{d1} and \code{d2} must be positive integer. 
#' For all other methods, the inputs must be a dummy variable indication 0=Nay and 1=Yea (and 9 for unit nonresponse).
#' All methods allow item non-response as NA.
#' @param dims Number of dimension in ideal point computation. Must be a positive 
#' integer between 1 and 10. If \code{bridge.method=="homography"}, only 2 is accepted 
#' at this point.  
#' @param polarity A vector specifying the row number of the \code{d1}(or pooled data)'s respondent(s) 
#' constrained to have a positive (i.e., right-wing or conservative) score on each dimension.
#' Used when \code{method} is \code{"ooc"}, \code{"oc"}, or \code{"wnominate"}. 
#' @param irt.type Used if \code{method} is \code{"irtMCMC"} or \code{"irtEM"}. 
#' If \code{"binary"} the binary IRT will be estimated. If \code{"ordered"}, 
#' the ordered IRT will be estimated (the ordered option is provisional). 
#' @param ... Additional arguments passed to the ideal point estimation function 
#' called in \code{method}.
#'   
#' @return A list with the following elements
#' \itemize{
#'   \item \code{ip}: Matrix of idealpoints.
#'   \item \code{info}: Outputs including other information.
#' }
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' 
#' @seealso \code{\link{ipbridging}}, \code{\link{oocflex}}, \code{\link[oc]{oc}}, 
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
#' 
#' @export

ipest <- function(d, 
                  method="ooc", 
                  dims=2,
                  polarity = c(1,1),
                  irt.type = "binary",
                  ...) {
  
  ## Check Inputs
  if (is.vector(d)) {
    warning("d is a vector. Converted to matrix with one column.")
    d <- as.matrix(d, ncol=1)
  }
  if (!is.matrix(d)&!is.data.frame(d)) stop("d must be matrix or data.frame!") 
  if (is.data.frame(d)) d <- as.matrix(d)
  
  # Estimate Ideal Points
  if (method=="ooc") {
    
    ip_f <- oocflex(d, dims=dims, polarity=polarity, ...)
    ip <- ip_f$respondents[,grepl("coord", colnames(ip_f$respondents))]
    
  } else if (method=="oc") {
    
    # Default minvotes set to 10
    set.minvotes <- TRUE
    if (length(list(...))>0) {
      if ("minvotes"%in%names(list(...))) set.minvotes <- FALSE
    }  
    
    if (set.minvotes) {
      ip_f <- oc(rollcall(d), dims=dims, polarity=polarity, 
                 minvotes = 10, ...)
    } else {
      ip_f <- oc(rollcall(d), dims=dims, polarity=polarity, ...)
    }
    ip <- ip_f$legislators[,grepl("coord", colnames(ip_f$legislators))]
    
  } else if (method=="wnominate") {
    
    # Default minvotes set to 10
    set.minvotes <- TRUE
    if (length(list(...))>0) {
      if ("minvotes"%in%names(list(...))) set.minvotes <- FALSE
    }  
    
    if (set.minvotes) {
      ip_f <- wnominate(rollcall(d), dims=dims, polarity=polarity, 
                        minvotes=10, ...)
    } else {
      ip_f <- wnominate(rollcall(d), dims=dims, polarity=polarity, ...)
    }
    ip <- ip_f$legislators[,grepl("coord", colnames(ip_f$legislators))]
    
  } else if (method=="irtMCMC") {
    
    if (irt.type=="binary") {
      ip_f <- ideal(rollcall(d), d=dims, ...)
      ip <- ip_f$xbar
    } else if (irt.type=="ordered") {
      ip_f <- MCMCpack::MCMCfactanal(d, factors=dims, ...)
      ip <- "NOT YET RETRIEVED"
    }
    
  } else if (method=="irtEM") {
    
    if (irt.type=="binary") {
      ip_f <- emIRT::binIRT(d, ...)
      ip <- ip_f$means$x
    } else if (irt.type=="ordered") {
      ip_f <- emIRT::ordIRT(d, ...)
      ip <- ip_f$means$x
    }
    
  } else if (method=="blackbox") {
    
    # Check if minscale argument is manually set
    # *minscale has no default in blackbox function
    set.minscale <- TRUE
    if (length(list(...))>0) {
      if ("minscale"%in%names(list(...))) set.minscale <- FALSE
    }  
    
    d <- as.data.frame(apply(d, 2, as.numeric))
    
    # Blackbox Scaling
    if (set.minscale) {
      # Set minscale to 10 if not manually set
      ip_f <- blackbox(d, dims=dims, minscale=10, ...)
    } else {
      ip_f <- blackbox(d, dims=dims, ...)
    }
    ip <- ip_f$individuals[[dims]]
    
  } else {
    
    stop("Invalid 'method' value!")
    
  }
  
  out <- list(ip=ip, model=ip_f)  
  # Additional Arguments
  out$method = method
  out$dims = dims
  out$polarity = polarity
  if (method%in%c("irtMCMC","irtEM")) out$irt.type=irt.type
  
  class(out) <- c("ipest.out", class(out))
  
  return(out)
  
}
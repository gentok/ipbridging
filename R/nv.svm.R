#' Estimating Normal Vector through Support Vector Machine
#' 
#' @description Compute normal vector from OC coordinates and 
#' ordered responses.
#' 
#' @param xmat Matrix of OC coordinates (i.e., predictors).
#' @param resp Response Variable (i.e., ordered choices).
#' @param kernel The kernel used in training and predicting. The 
#' default is \code{"radial"} to use radial basis function (RBF). 
#' The alterantaives are \code{"linear"}, \code{"polynomial"} or \code{"sigmoid"}. 
#' See the help file of \code{\link[e1071]{svm}} for more details.
#' @param tune.param Method to determine the parameters. 
#' Following options are currently available:
#' \itemize{
#'   \item \code{NULL}: Use \code{\link[e1071]{svm}} with the default 
#'   or manually set parameter values. To manually set parameters.
#'   Check the available parameters (and their default values) in 
#'   \code{\link[e1071]{svm}} function in \code{e1071} package.
#'   \item \code{"heuristics"} (default): If \code{kernel=="radial"}, 
#'   Use heuristics-based method (i.e., formulation identical to 
#'   \code{\link[kernlab]{sigest}} in \code{kernlab} package) 
#'   to determine optimal gamma. Additionally, if \code{resp} is numeric, 
#'   use heuristic based method described in Cherkassky & Ma (2004) to 
#'   determine optimal parameters for C (cost) and epsilon (epsilon-regression is 
#'   automatically selected).
#'   If \code{resp} is a factor, the same behavior as \code{NULL}. 
#'   \item A named list of parameter vectors spanning the sampling space. 
#'   If list is assigned here, the method use \code{\link[e1071]{best.tune}} function
#'   to obtain the best perfomed model in the given parameter ranges. 
#' }
#' @param param.heuristics.k If \code{tune.param=="heuristics"}, the method uses 
#' k nearest neighbors regression to estimate the noise variance in response variable
#'  for the detemination of epsilon.
#' This argument sets k for this method. The default is 3. Cherkassky & Ma (2004) 
#' recommends somewhere between 2 to 6. 
#' @param param.heuristics.frac Fraction of data to use for heuristics-based estimation of gamma. 
#' By default, 50 percent of the data is used to estimate the range of the gamma hyperparameter.
#' @param ... Additional arguments passed to \code{\link[e1071]{svm}} 
#' or \code{\link[e1071]{best.tune}}.
#' 
#' @return A vector of coefficients.
#' 
#' @seealso \code{\link[e1071]{svm}}, \code{\link[e1071]{best.tune}}
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' This is a modified version of the code included in \code{\link[ooc]{ooc}} function.
#' 
#' @references 
#' \itemize{
#'   \item Chang, C. & C. Lin, LIBSVM: A Library for Support Vector Machines, \url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}.
#'   \item Cherkassky, V. & Y. Ma, 2004, "Practical Selection of SVM Parameters and Noise Estimation for SVM Regression", Neural Networks, 17, 113 - 126.
#' }
#' 
#' @export

nv.svm <- function(xmat, resp, 
                   kernel="radial",
                   tune.param = "heuristics", 
                   param.heuristics.k = 3, 
                   param.heuristics.frac = 0.5, ...) {
  
  # cat(paste0("\nxmat dimension is ",dim(xmat)[1]," rows and ", dim(xmat)[2], " columns. ", "resp length is ", length(resp), ".\n\n"))
  
  if (is.null(tune.param)) {
    # No Tuning
    
    res <- e1071::svm(x=xmat, y=resp, scale=FALSE, kernel=kernel, ...)
    coefs <- t(res$coefs) %*% res$SV
    
  } else if (is.list(tune.param)) {
    # Tuning Parameter Values Specified  
    
    res <- e1071::best.tune(e1071::svm, train.x=xmat, train.y=resp, 
                            ranges = tune.param, scale=FALSE, 
                            kernel=kernel, ...)
    coefs <- t(res$coefs) %*% res$SV
    
  } else if (tune.param=="heuristics") {
    
    # Drop missing cases
    respcp <- resp[complete.cases(xmat & resp)]
    xmatcp <- xmat[complete.cases(xmat & resp),,drop=FALSE]
    
    # Optimal Gamma
    if (kernel=="radial") {
      # Use heuristics to set gamma (follows sigest() function from kernlab package)
      # Refer to sigest() function in kernlab package
      m <- dim(xmatcp)[1]
      n <- floor(param.heuristics.frac*m)
      index <- sample(1:m, n, replace = TRUE)
      index2 <- sample(1:m, n, replace = TRUE)
      temp <- xmatcp[index,, drop=FALSE] - xmatcp[index2,,drop=FALSE]
      dist <- rowSums(temp^2)
      optg <- 1/quantile(dist[dist!=0],probs=c(0.5)) 
      # Can be any value between 10-90%tile but choose median
    } else {
      optg <- 1/ncol(xmatcp) # Default for e1071::svm
    }
    
    if (is.numeric(resp)) {
      # Use Heuristics to set epsilon and C if using SVM EPS regression  
      # Use eps regression and set C and epsilon by heuristics 
      # described in Cherkassky & Ma (2004)
      
      # Optimal C
      optc <- max( abs(c(mean(respcp)+3*sd(respcp))),  abs(c(mean(respcp)-3*sd(respcp))) )
      # Optimal epsilon
      resphat <- predict(caret::knnreg(x=xmatcp, y=respcp, k=param.heuristics.k), xmatcp)
      tmp <- (length(respcp)^(1/5))*param.heuristics.k
      s <- sqrt( (tmp/(tmp-1)) * (1/length(respcp)) * sum((respcp-resphat)^2))
      opteps <- 3*s*sqrt(log(length(respcp))/length(respcp))
      # Run EPS regression with Optimal C and Epsilon
      res <- e1071::svm(x=xmat, y=resp, scale=FALSE, 
                        type="eps-regression", cost=optc, epsilon=opteps, 
                        kernel=kernel, gamma=optg, ...)
      coefs <- t(res$coefs) %*% res$SV
      
    } else if (is.factor(resp)) {
      # Only set gamma if resp is a factor
      
      res <- e1071::svm(x=xmat, y=resp, scale=FALSE, 
                        kernel=kernel, 
                        gamma=optg, ...)
      coefs <- t(res$coefs) %*% res$SV
      
    } else {
      stop("Invalid 'resp' class!")
    }
    
  } else {
    stop("Invalid 'tune.param' argument!")
  }
  
  return(coefs)
  
}

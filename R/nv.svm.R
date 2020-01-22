#' Estimating Normal Vector through Support Vector Machine
#' 
#' @param xmat Matrix of OC coordinates (i.e., predictors).
#' @param resp Response Variable (i.e., ordered choices).
#' @param tune.param Method to determine the parameters. 
#' Following options are currently available:
#' \itemize{
#'   \item \code{NULL}: Use \code{\link[e1071]{svm}} with the default 
#'   or manually set parameter values. To manually set parameters.
#'   Check the available parameters (and their default values) in 
#'   \code{\link[e1071]{svm}} function in \code{e1071} package.
#'   \item \code{"heuristics"} (default): If \code{resp} is numeric, 
#'   use heuristic based method described in Cherkassky & Ma (2004) to 
#'   determine optimal parameters for C (cost) and epsilon.
#'   If \code{resp} is a factor, the same behavior as \code{NULL}. 
#'   \item A named list of parameter vectors spanning the sampling space. 
#'   If list is assigned here, the method use \code{\link[e1071]{best.tune}} function
#'   to obtain the best perfomed model in the given parameter ranges. 
#' }
#' @param use.gpu If \code{TRUE}, the method replaces \code{\link[e1071]{best.tune}} 
#' and \code{\link[e1071]{svm}} functions with 
#' the \code{best.tune} and \code{svm} functions in \code{Rgtsvm} package, which contains
#' GPU accelerate method of \code{svm} function. Install \code{Rgtsvm} manually. 
#' You need to have Linux computer with NVIDIA GPU for the installation.
#' Check \url{https://github.com/Danko-Lab/Rgtsvm} for more details.  
#' @param param.heuristics.k If \code{tune.param=="heuristics"}, the method uses 
#' k nearest neighbors regression to estimate the noise variance in response variable.
#' This argument sets k for this method. The default is 3. Cherkassky & Ma (2004) 
#' recommends somewhere between 2 to 6. 
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
#'   \item Cotter, A., N. Srebro, & J. Kreshet, 2011, "A GPU-Tailored Approach for Training Kernelized SVMs", 17th ACM SIGKDD Conference on Knowledge Discovery and Data Mining.
#' }
#' 
#' @export

nv.svm <- function(xmat, resp, 
                   tune.param = "heuristics", 
                   use.gpu = FALSE, 
                   param.heuristics.k = 3, ...) {
  
  # cat(paste0("\nxmat dimension is ",dim(xmat)[1]," rows and ", dim(xmat)[2], " columns. ", "resp length is ", length(resp), ".\n\n"))
  
  if (use.gpu==FALSE) {
  # NOT using GPU
      
    if (is.null(tune.param)) {
      # No Tuning
      
      res <- e1071::svm(x=xmat, y=resp, scale=FALSE, ...)
      coefs <- t(res$coefs) %*% res$SV
      
    } else if (is.list(tune.param)) {
      # Tuning Parameter Values Specified  
      
      res <- e1071::best.tune(e1071::svm, train.x=xmat, train.y=resp, 
                              ranges = tune.param, scale=FALSE, ...)
      coefs <- t(res$coefs) %*% res$SV
      
    } else if (tune.param=="heuristics") {
      # Use Heuristics to set epsilon and C if using SVM EPS regression  
      
      if (is.numeric(resp)) {
        # Use eps regression and set C and epsilon by heuristics 
        # described in Cherkassky & Ma (2004)
        
        # Drop missing cases
        respcp <- resp[complete.cases(xmat & resp)]
        xmatcp <- xmat[complete.cases(xmat & resp),]
        # Optimal C
        optc <- max( abs(c(mean(respcp)+3*sd(respcp))),  abs(c(mean(respcp)-3*sd(respcp))) )
        # Optimal epsilon
        resphat <- predict(caret::knnreg(x=xmatcp, y=respcp, k=param.heuristics.k), xmatcp)
        tmp <- (length(respcp)^(1/5))*param.heuristics.k
        s <- sqrt( (tmp/(tmp-1)) * (1/length(respcp)) * sum((respcp-resphat)^2))
        opteps <- 3*s*sqrt(log(length(respcp))/length(respcp))
        # Run EPS regression with Optimal C and Epsilon
        res <- e1071::svm(x=xmat, y=resp, scale=FALSE, 
                          type="eps-regression", cost=optc, epsilon=opteps, ...)
        coefs <- t(res$coefs) %*% res$SV
        
      } else if (is.factor(resp)) {
        # Use Default (or manually set) parameters if resp is a factor
        
        res <- e1071::svm(x=xmat, y=resp, scale=FALSE, ...)
        coefs <- t(res$coefs) %*% res$SV
        
      } else {
        stop("Invalid 'resp' class!")
      }
      
    } else {
      stop("Invalid 'tune.param' argument!")
    }

  } else {
  ## USE GPU
      
    if (!requireNamespace("Rgtsvm", quietly = TRUE)) {
      stop("Package \"Rgtsvm\" needed for this option to work. Please install it.",
           call. = FALSE)
    }
    
    settype <- ifelse(is.numeric(resp),"eps-regression", "C-classification")
    
    if (is.null(tune.param)) {
      # No Tuning
      
      res <- Rgtsvm::svm(x=xmat, y=resp, scale=FALSE, type=settype, ...)
      coefs <- t(res$coefs[seq(1,nrow(res$SV))]) %*% res$SV
      
    } else if (is.list(tune.param)) {
      # Tuning Parameter Values Specified  
      
      ## SVM with Tuning, using GPU through Rgtsvm package (only Linux)
      
      res <- Rgtsvm::tune.svm(Rgtsvm::svm, x=xmat, y=resp, 
                               degree=tune.param$degree,
                               gamma=tune.param$gamma,
                               coef0=tune.param$coef0,
                               cost=tune.param$cost,
                               nu=tune.param$nu,
                               class.weights=tune.param$class.weights,
                               epsilon=tune.param$epsilon, 
                               type = settype, scale=FALSE, ...)
      res <- res$best.model
      coefs <- t(res$coefs[seq(1,nrow(res$SV))]) %*% res$SV
      
    } else if (tune.param=="heuristics") {
      # Use Heuristics to set epsilon and C if using SVM EPS regression  
      
      if (is.numeric(resp)) {
        # Use eps regression and set C and epsilon by heuristics 
        # described in Cherkassky & Ma (2004)
        
        # Drop missing cases
        respcp <- resp[complete.cases(xmat & resp)]
        xmatcp <- xmat[complete.cases(xmat & resp),]
        # Optimal C
        optc <- max( abs(c(mean(respcp)+3*sd(respcp))),  abs(c(mean(respcp)-3*sd(respcp))) )
        # Optimal epsilon
        resphat <- predict(caret::knnreg(x=xmatcp, y=respcp, k=param.heuristics.k), xmatcp)
        tmp <- (length(respcp)^(1/5))*param.heuristics.k
        s <- sqrt( (tmp/(tmp-1)) * (1/length(respcp)) * sum((respcp-resphat)^2))
        opteps <- 3*s*sqrt(log(length(respcp))/length(respcp))
        # Run EPS regression with Optimal C and Epsilon
        res <- Rgtsvm::svm(x=xmat, y=resp, scale=FALSE, 
                          type="eps-regression", cost=optc, epsilon=opteps, ...)
        coefs <- t(res$coefs[seq(1,nrow(res$SV))]) %*% res$SV
        
      } else if (is.factor(resp)) {
        # Use Default (or manually set) parameters if resp is a factor
        
        res <- Rgtsvm::svm(x=xmat, y=resp, scale=FALSE, type=settype, ...)
        coefs <- t(res$coefs[seq(1,nrow(res$SV))]) %*% res$SV
        
      } else {
        stop("Invalid 'resp' class!")
      }
      
    } else {
      stop("Invalid 'tune.param' argument!")
    }
    
  }

  return(coefs)
  
}

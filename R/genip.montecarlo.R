#' Generate Simulated Responses of Respondents with Ideal Points
#' 
#' @description Generate simulated responses of survey respondents with ideal points 
#' through Monte Carlo method.
#' 
#' @param n Number of respondents.
#' @param q Number of issues.
#' @param ncat Number of ordered category in responses.
#' @param ndim Number of dimentions in ideal points.
#' @param missing Proportion of missing responses. Ranges from 0-1.
#' @param g1.size Size of the first group defined as proportion (0-1 range). 
#' The second group's size is \code{1-g1.size}. 
#' @param correlations.lim.g1 Limits in the range of correlation between respondents and issues in the group 1. 
#' Higher value indicates stronger association between respondent's values and their 
#' ideal points in the issue.
#' @param correlations.lim.g2 Limits in the range of correlation between respondents and issues in the group 2. 
#' @param utility.diff.level Must be a positive integer, default is 1. 
#' In this simulation, the value of \code{correlations} is potentially related to the 
#' likely set of utility function form, that converts ideal points into ordered survey 
#' responses. \code{utility.diff.level} deterimines the extent of differences in likely 
#' utility function forms (randomly selected from linear, normal, and qudratic utility functions, 
#' with some probability weights). The value of 1 indicates no relationship between 
#'\code{correlations} and likely utility function forms, and higher values indicate 
#' more extreme relationship between \code{correlations} and likely utility function forms.   
#' @param error.respondents The lower and upper limits of uniform distribution that are used to initiate values of respondents.
#' @param error.issues The \code{shape} and \code{rate} parameter in \code{\link[stats]{rgamma}} distribution, respectively. Used to set initial values of issue positions. 
#' @param idealpoints.lim Limits in the range of ideal points. The values outside 
#' of this range is recoded into minimum and maximum value of the range.
#' 
#' @return A list with the following elements
#' \itemize{
#'   \item \code{simulated.responses}: n*q data.frame of simulated responses.
#'   \item \code{perfect.responses}: Hypothetical responses if each respondents give responses without error.
#'   \item \code{idealpoints}: Ideal point coordinates of each respondent.
#'   \item \code{normalvectors}: Normal vectors.
#'   \item \code{heteroskedastic.respondents}: Initial values assigned to respondents.
#'   \item \code{heteroskedastic.issues}: Initial values assigned to issues.
#'   \item \code{correlations}: Correlation between ideal point dimensions.
#'   \item \code{knowledge}: Correlation binned into three categories.
#'   \item \code{groups}: Group identifiers.
#'   \item \code{error}: Proportions of incorrect choices.
#' }
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' 
#' @import stats
#' 
#' @export

genip.montecarlo <- function(n=1200, 
                           q=20,
                           ncat=5,
                           ndim=2,
                           missing=0.1,
                           g1.size = 0.5,
                           utility.diff.level = 1,
                           correlations.lim.g1=c(-0.1,0.7), 
                           correlations.lim.g2=c(-0.1,0.7), 
                           error.respondents=c(0.1,1), 
                           error.issues=c(2,1),
                           idealpoints.lim=c(-2,2)){
  
  # 0.) Generate Respondents and Issues
  ## Respondents
  heteroskedastic.respondents <- runif(n, error.respondents[1], error.respondents[2])
  ## Issues
  heteroskedastic.issues <- rgamma(q, error.issues[1], error.issues[2])
  ## Groups
  groups <- c(rep(1,floor(n*g1.size)),rep(2,n-floor(n*g1.size))) 
  ## Correlations b/w Each Respondent and Issues
  correlations.g1 <- runif(floor(n*g1.size), correlations.lim.g1[1], correlations.lim.g1[2])
  correlations.g2 <- runif(n - floor(n*g1.size), correlations.lim.g2[1], correlations.lim.g2[2])
  correlations <- c(correlations.g1, correlations.g2)
  knowledge <- statar::xtile(correlations, 3)
  #
  # 1.) Generate respondent ideal points
  #
  mu <- rep(0, ndim)
  Sigma <- list()
  idealpoints <- matrix(NA, nrow=n, ncol=ndim)
  for (i in 1:n){
    Sigma[[i]] <- matrix(1, nrow=ndim, ncol=ndim)
    Sigma[[i]][lower.tri(Sigma[[i]])] <- runif(sum(lower.tri(Sigma[[i]])), correlations[i]-0.2, correlations[i]+0.2)
    Sigma[[i]][upper.tri(Sigma[[i]])] <- t(Sigma[[i]])[upper.tri(t(Sigma[[i]]))]
    Sigma[[i]] <- as.matrix(Matrix::nearPD(Sigma[[i]])$mat)
    idealpoints[i,] <- MASS::mvrnorm(1, mu=mu, Sigma=Sigma[[i]])
  }
  #
  idealpoints[idealpoints > idealpoints.lim[2]] <- idealpoints.lim[2]
  idealpoints[idealpoints < idealpoints.lim[1]] <- idealpoints.lim[1]
  #
  # 2.) Generate normal vectors
  #
  normalvectors <- hitandrun::hypersphere.sample(ndim, q)
  #
  for (j in 1:q){
    if (normalvectors[j,1] < 0) normalvectors[j,] <- -1 * normalvectors[j,]
  }
  #
  # 3.) Project respondents on normal vectors
  #
  respondent.projections <- idealpoints %*% t(normalvectors)
  #
  # 4.) Generate outcome locations
  #
  outcome.locations <- apply(respondent.projections, 2, function(x){sort(runif(ncat,min(x),max(x)))})
  #
  # 5.) Define utility functions
  #
  linear.utility <- function(idealpoint, choices){
    tmp <- 1 - (3*abs(idealpoint - choices))
    return(tmp)}
  #
  normal.utility <- function(idealpoint, choices){
    tmp <- 5 * exp(-1 * ((idealpoint - choices)^2))
    return(tmp)}
  #
  quad.utility <- function(idealpoint, choices){
    tmp <- 1 - 2 * (idealpoint - choices)^2
    return(tmp)}
  #
  # 6.) Calculate utility and response probabilities
  #
  ## Sampling two sets of utility weights 
  utmp <- sample.int(utility.diff.level,6,replace=TRUE)
  # Correlations in 0-1 Range
  rtmp <- range(c(correlations.lim.g1, correlations.lim.g2))
  stdcors <- (correlations - rtmp[1])/(rtmp[2]-rtmp[1])
  # Correlations and Utility Function is Related
  uprobs <- sapply(1:3, function(k) utmp[k]*stdcors + utmp[k+3]*(1-stdcors))
  # 
  utility.fun <- 
    apply(uprobs, 1, function(k) sample(c("linear","normal","quadratic"), 1, prob=k))
  # utility.fun <- sample(c("linear","normal","quadratic"), n, 
  #                       replace=TRUE, prob=utility.probs)
  #
  systematic.utility <- total.utility <- probmat <- list()
  for (j in 1:q){
    systematic.utility[[j]] <- matrix(NA, nrow=n, ncol=ncat)
    total.utility[[j]] <- matrix(NA, nrow=n, ncol=ncat)
    probmat[[j]] <- matrix(NA, nrow=n, ncol=ncat)
  }
  #
  for (i in 1:n){
    for (j in 1:q){
      #
      if(utility.fun[i]=="linear"){
        systematic.utility[[j]][i,] <- linear.utility(respondent.projections[i,j], outcome.locations[,j])
      }
      #
      if(utility.fun[i]=="normal"){
        systematic.utility[[j]][i,] <- normal.utility(respondent.projections[i,j], outcome.locations[,j])
      }
      #
      if(utility.fun[i]=="quadratic"){
        systematic.utility[[j]][i,] <- quad.utility(respondent.projections[i,j], outcome.locations[,j])
      }
      #
      total.utility[[j]][i,] <- systematic.utility[[j]][i,] / exp(heteroskedastic.respondents[i] * heteroskedastic.issues[j])
      probmat[[j]][i,] <- pnorm(total.utility[[j]][i,]) / sum(pnorm(total.utility[[j]][i,]))
    }}
  #
  for (j in 1:q){
    probmat[[j]][is.na(probmat[[j]])] <- 1
  }
  #
  # 7.) Generate perfect and error (probabilistic) responses
  #
  simulated.responses <- perfect.responses <- matrix(NA, nrow=n, ncol=q)
  #
  for (i in 1:n){
    for (j in 1:q){
      simulated.responses[i,j] <- sample(1:ncat, 1, prob=probmat[[j]][i,])
      # Note: perfect voting could of course also be simulated using the maximum of total.utility
      perfect.responses[i,j] <- which.max(systematic.utility[[j]][i,])
    }}
  #
  # 8.) Insert missing data at random
  #
  if(missing > 0){
    miss <- expand.grid(n=1:n, q=1:q)
    missing.selected <- as.matrix(miss)[sample(1:nrow(miss), (n*q*missing), replace=FALSE),]
    for (i in 1:nrow(missing.selected)){
      simulated.responses[missing.selected[i,1], missing.selected[i,2]] <- NA
      perfect.responses[missing.selected[i,1], missing.selected[i,2]] <- NA
    }
  }
  #
  # 9.) Organize output
  #
  correctvotes <- sum(diag(table(simulated.responses,perfect.responses)))
  totalvotes <- sum(table(simulated.responses,perfect.responses))
  #
  return(list(simulated.responses = simulated.responses,
              perfect.responses = perfect.responses,
              idealpoints = idealpoints,
              normalvectors = normalvectors,
              heteroskedastic.respondents = heteroskedastic.respondents,
              heteroskedastic.issues = heteroskedastic.issues,
              correlations = correlations,
              knowledge = knowledge,
              groups = groups,
              error = (1 - (correctvotes / totalvotes))
  ))
  #
}

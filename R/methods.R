#' Ordered Optimal Classification Summary
#' 
#' @description \code{summary.oocflexObject} reads an Ordered Optimal 
#' Classification object exported by \code{\link{oocflex}} function and
#' prints a summary.
#' 
#' @param object an \code{oocflexObject} output object.
#' @param verbose logical, includes all ideal points if TRUE, 
#' otherwise only returns the first 10 respondents.
#' @param ... Other arguments. Do nothing and are not passed to any functions.
#' 
#' @return A summary of a \code{oocflexObject}. 
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' This code is the modified version of the \code{\link[oc]{summary.OCobject}} 
#' method written by Keith Poole, Jeffrey Lewis, James Lo, and Royce Carroll.
#' 
#' @seealso \code{\link{oocflex}}
#' 
#' @export

summary.oocflexObject <- function (object, verbose = FALSE, ...) 
{
  if (!class(object) == "oocflexObject") 
    stop("Input is not of class 'oocflexObject'.")
  cat("\n\nSUMMARY OF ORDERED OPTIMAL CLASSIFICATION")
  cat("\n----------------------------\n")
  cat("\nNumber of Respondents:\t ", dim(na.omit(object$respondents))[1], 
      " (", dim(object$respondents)[1] - dim(na.omit(object$respondents))[1], 
      " respondents deleted)", sep = "")
  cat("\nNumber of Issues:\t ", dim(na.omit(object$issues))[1], 
      " (", dim(object$issues)[1] - dim(na.omit(object$issues))[1], 
      " issues deleted)", sep = "")
  cat("\nNumber of Dimensions:\t", object$dimensions)
  correctYea <- sum(as.numeric(object$respondents[, "correctYea"]), 
                    na.rm = TRUE)
  allYea <- correctYea + sum(as.numeric(object$respondents[, 
                                                           "wrongNay"]), na.rm = TRUE)
  correctNay <- sum(as.numeric(object$respondents[, "correctNay"]), 
                    na.rm = TRUE)
  allNay <- correctNay + sum(as.numeric(object$respondents[, 
                                                           "wrongYea"]), na.rm = TRUE)
  cat("\nPredicted Yeas:\t\t ", correctYea, " of ", allYea, 
      " (", round(100 * correctYea/allYea, 1), "%) predictions correct", 
      sep = "")
  cat("\nPredicted Nays:\t\t ", correctNay, " of ", allNay, 
      " (", round(100 * correctNay/allNay, 1), "%) predictions correct\n\n", 
      sep = "")
  if (!verbose) {
    cat("The first 10 respondent estimates are:\n")
    if (object$dimensions != 1) {
      round(object$respondents[1:10, paste("coord", 1:object$dimensions, 
                                           "D", sep = "")], 3)
    }
    else {
      round(object$respondents[1:10, c("coord1D"), drop = FALSE], 
            3)
    }
  }
  else {
    if (object$dimensions != 1) {
      round(object$respondents[, paste("coord", 1:object$dimensions, 
                                       "D", sep = "")], 3)
    }
    else {
      round(object$respondents[, c("coord1D"), drop = FALSE], 
            3)
    }
  }
}
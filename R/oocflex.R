#' Compute Normal Vector
#' 
#' @param xmat Matrix of OC coordinates (i.e., predictors).
#' @param resp Response Variable (i.e., ordered choices).
#' @param nv.method Method of normal vector estimation. 
#' Following pre-defined methods are currently available:
#' \itemize{
#'   \item \code{"oprobit"}: Ordered probit regression (calls \code{\link[MASS]{polr}} function from \code{MASS} package)
#'   \item \code{"ologit"}: Ordered logistic regression (calls \code{\link[MASS]{polr}} function from \code{MASS} package)
#'   \item \code{"svm.reg"}: Support Vector Regression (calls \code{\link{nv.svm}} function)
#'   \item \code{"svm.class"}: Support Vector Classification with dichotomized response (calls \code{\link{nv.svm}} function)
#'   \item \code{"krls"}: Kernel Regularized Least Squares (calls \code{\link{nv.krls}} function). 
#' }
#' Alternatively, one can set custom made function that (at least) feeds in 
#' \code{xmat} and \code{resp} arguments and returns the vector of coefficients 
#' (the length must be equal to the number of columns of \code{xmat}). 
#' 
#' @param ... Additional arguments passed to the function assigned by \code{nv.method}.
#' 
#' @return The vector of normal vector estimates.
#' 
#' @author Tzu-Ping Liu \email{jamesliu0222@@gmail.com}, Gento Kato \email{gento.badger@@gmail.com}, and Sam Fuller \email{sjfuller@@ucdavis.edu}.
#' 
#' @export

compute.nv <- function(xmat, resp, nv.method, ...) {
  
  if ("function"%in%class(nv.method)) {
    
    coefs <- nv.method(xmat=xmat, resp=resp, ...)
    
  } else {
    
    if (!is.character(nv.method)|length(nv.method)!=1) {
      
      stop("Invalid value assigned to nv.method!")
      
    } else if (nv.method=="oprobit"){
      
      ## Ordinal Probit
      coefs <- MASS::polr(factor(resp) ~ xmat, method="probit", ...)$coefficients
      
    } else if (nv.method=="ologit"){
      
      ## Ordinal Logit
      coefs <- MASS::polr(factor(resp) ~ xmat, method="probit", ...)$coefficients
      
    } else if (nv.method%in%c("svm.reg","svm.class")){
      
      ## SVM Regression, numeric response
      if (nv.method=="svm.reg") resp <- as.numeric(resp)
      ## SVM Classification, dichotomize response
      if (nv.method=="svm.class") resp <- factor(ooc::binary.balanced(resp))
      ## Get Normal Vector through SVM
      coefs <- nv.svm(xmat=xmat, resp=resp, ...) 

    } else if (nv.method=="krls"){
      
    # Kernel Regularized Least Squares
    coefs <- nv.krls(xmat=xmat, resp=resp, ...) 
    
    }
    
  }
  
  ## Normal Vector Estiamtes  
  est <- coefs / sqrt(sum(coefs^2))
  
  return(est)
}

#' Ordered Optimal Classification (Flexible Version)
#' 
#' @description Performs Ordered Optimal Classification (OOC) with flexible strategies of estimation. OOC is an extension of Poole's (2000) nonparametric unfolding procedure for the analysis of ordinal choice data (e.g., “Strongly Agree,” “Somewhat Agree,” “Somewhat Disagree,” “Strongly Disagree”).
#' 
#' @param votemat A matrix of ordinal choice data, must be consecutive 
#' integers starting with 1. Binary choice should be coded so that 
#' 1 = Support, 2 = Oppose. Missing data must be coded as NA.
#' @param dims Number of dimensions to estimate.
#' @param minvotes Minimum number of votes required for a respondent to be 
#' included in the analysis.
#' @param lop A proportion between 0 and 1, the cut-off used for excluding 
#' lopsided votes.
#' @param polarity A vector specifying the row number of the respondent(s) 
#' constrained to have a positive (i.e., right-wing or conservative) score 
#' on each dimension.
#' @param iter Number of iterations of the modified Optimal Classification 
#' algorithm.
#' @param nv.method The method used to compute normal vectors at each step 
#' of the iteration. Current choices include: 
#' \code{"oprobit"} (Ordered probit regression), 
#' \code{"ologit"} (Ordered logistic regression), 
#' \code{"svm.reg"} (SVM: regression), 
#' \code{"svm.class"} (SVM: classification), 
#' and \code{"krls"} (Kernel regularized least squares). 
#' One can also assign a custom-made function to compute normal vectors.
#' See more details in \code{\link{compute.nv}}.
#' @param binary.method The method used to transform ordinal responses 
#' to binary ones. Following options are currently available:
#' \itemize{
#'   \item \code{"dominance"} (default): Transformed to the matrix with binary choices organized in a dominance pattern.
#'   \item \code{"even"}: Transformed to the matrix with binary choices recoded to be as evenly-balanced as possible.
#' }
#' @param random.seed If not \code{NULL}, scalar indicates the random seed for reproducibility.
#' @param ... Additional arguments passed to the function assigned 
#' by \code{nv.method} (See more details in \code{\link{compute.nv}}). 
#' 
#' @return An \code{oocflexObject} with the following elements
#' \itemize{
#'   \item \code{respondents}: A matrix containing the respondent estimates.
#'   \item \code{issues}: A matrix containing the estimates for the (c-1) categories of each issue.
#'   \item \code{issues.unique}: A matrix containing only the estimates for the unique issues.
#'   \item \code{fits}: The percentage of all binary choices correctly predicted, 
#'   the Aggregate Proportional Reduction in Error (APRE) statistic for the binary choices, 
#'   the percentage of all ordinal choices correctly predicted, 
#'   and the Aggregate Proportional Reduction in Error (APRE) statistic for 
#'   the ordinal choices.
#'   \item \code{OC.result.binary}: Standard (binary) OC results 
#'   for the choice matrix arranged in a binary dominance pattern.
#'   \item \code{errorcount}: Error counts at each iteration of the Optimal Classification algorithm.
#'   \item \code{votemat.binary}: The matrix of ordinal choices recoded into binary categories 
#'   that are as balanced as possible.
#' }
#' 
#' @seealso \code{\link{compute.nv}}, \code{\link{nv.svm}}, and \code{\link{nv.krls}}
#' 
#' @author Tzu-ping Liu \email{jamesliu0222@@gmail.com} and Gento Kato \email{gento.badger@@gmail.com}. 
#' This code is the modified version of the \code{\link[ooc]{ooc}} 
#' function written by Christopher Hare and Keith T. Poole.
#' 
#' @references 
#' \itemize{
#'   \item Hare, C., Liu, T., & Lupton, R. N. 2018. "What Ordered Optimal Classification Reveals about Ideological Structure, Cleavages, and Polarization in the American Mass Public", Public Choice, 176, 57-78. 
#'   \item Armstrong, David A., Ryan Bakker, Royce Carroll, Christopher Hare, Keith T. Poole, and Howard Rosenthal. 2014. Analyzing Spatial Models of Choice and Judgment with R. Boca Raton, FL: CRC Press.
#'   \item Poole, Keith T. 2000. "Nonparametric unfolding of BInary Choice Data." Political Analysis 8(e): 211-237.       
#' }
#' 
#' @import stats
#' @import oc
#' @import ooc
#' 
#' @export

oocflex <- function(votemat, 
                    dims=2, 
                    minvotes=10, 
                    lop=0.001, 
                    polarity=c(1,1),
                    iter=25, 
                    nv.method="svm.reg",
                    binary.method="dominance",
                    random.seed = NULL, ...){
  
  if (!is.null(random.seed)) set.seed(random.seed)
  # set.seed(1985)
  
  nvotes <- ncol(votemat) # Number of Issues
  ndim <- dims # Number of dimention in ideal points

  
  #  ***********************************************************************************************************************
  #  NOTE THAT ncol(votemat) is the number of input columns == number of issue scales being analyzed.  Each scale can have
  #    a different number of categories -- 3 Point Scales, 7 Point Scales, etc.  Hence, ZVEC[,s] is a set of stacked normal
  #    vectors that are the same for each issue. If it is a 3 Point Scale there will be two normal vectors; if it is a 7
  #    Point scale there will be six normal vectors; etc.
  #    Hence the length of ZVEC[,] = SUM_j=1,q{ncat_j-1}
  #
  #  **********************************************************************************************************************
  
  ## Number of choice categories for each question
  ncat <- NULL
  for (j in 1:nvotes){
    ncat[j] <- length(table(factor(votemat[,j], levels = min(votemat[,j], na.rm=TRUE):max(votemat[,j], na.rm=TRUE))))
  }
  sum.ncat <- sum(ncat)
  
  if (binary.method=="dominance") {
    
    # the matrix with binary choices organized in a dominance pattern
    votemat.binary <- do.call("cbind", plyr::alply(votemat, 2, binarize))

    #############################################################
    ## USE OPTIMAL CLASSIFICATION ON REGULAR MATRIX TO GET Nj  ##
    #############################################################
    hr.binary <- pscl::rollcall(votemat.binary, yea=1, nay=6, missing=9, notInLegis=NA)
    result.binary <- oc(hr.binary, dims=ndim, minvotes=minvotes, lop=lop, polarity=polarity)
    
    ################################################################
    ## THROW OUT RESPONDENTS AND ISSUES THAT DON'T MEET THRESHOLD ##
    ################################################################
    use.respondents <- !is.na(result.binary$legislators[,"coord1D"])
    votemat.binary <- votemat.binary[use.respondents,]
    votemat <- votemat[use.respondents,]
    nresp <- nrow(votemat.binary)
    
    ###############################################################  
    ## Running Ordered Optimal Classification (Single Dimension) ##
    ###############################################################
    if (ndim==1){
      
      respondents <- result.binary$legislators
      
      issues <- result.binary$rollcalls
      
      fits <- result.binary$fits
      
      minorityvote <- apply(cbind(issues[,1] + issues[,3], issues[,2] + issues[,4]),1,min)
      classerrors <- issues[,2] + issues[,3]
      correctclass <- issues[,1] + issues[,4]
      
      uniqueCLASSPERC <- sum(correctclass[1:(ncat[1]-1)]) / sum(correctclass[1:(ncat[1]-1)] + classerrors[1:(ncat[1]-1)])
      for (j in 2:nvotes){
        uniqueCLASSPERC[j] <- sum(correctclass[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]) / sum(correctclass[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])] + classerrors[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])])
      }
      
      uniquePRE <- sum(minorityvote[1:(ncat[1]-1)] - classerrors[1:(ncat[1]-1)]) / sum(minorityvote[1:(ncat[1]-1)])
      for (j in 2:nvotes){
        uniquePRE[j] <- sum(minorityvote[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])] - classerrors[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]) / sum(minorityvote[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])])
      }
      
      issues.unique <- cbind(uniqueCLASSPERC,uniquePRE)
      colnames(issues.unique) <- c("correctPerc","PRE")
      
      nvotes.binary <- nrow(issues)
      oc1 <- na.omit(respondents[,7])
      ws <- issues[,7]
      
      onedim.classification.errors <- matrix(NA,nrow=nresp,ncol=nvotes.binary)
      for (j in 1:nvotes.binary){
        ivote <- votemat.binary[,j]
        polarity <- oc1 - ws[j]
        errors1 <- ivote==1 & polarity >= 0
        errors2 <- ivote==6 & polarity <= 0
        errors3 <- ivote==1 & polarity <= 0
        errors4 <- ivote==6 & polarity >= 0
        kerrors1 <- ifelse(is.na(errors1),9,errors1)
        kerrors2 <- ifelse(is.na(errors2),9,errors2)
        kerrors3 <- ifelse(is.na(errors3),9,errors3)
        kerrors4 <- ifelse(is.na(errors4),9,errors4)
        kerrors12 <- sum(kerrors1==1)+sum(kerrors2==1)
        kerrors34 <- sum(kerrors3==1)+sum(kerrors4==1)
        if (kerrors12 >= kerrors34){
          yeaerror <- errors3
          nayerror <- errors4
        }
        if (kerrors12 < kerrors34){
          yeaerror <- errors1
          nayerror <- errors2
        }
        junk <- yeaerror
        junk[nayerror==1] <- 1
        onedim.classification.errors[,j] <- junk
      }
      
      correct.class.binary <- matrix(NA,nrow=nresp,ncol=nvotes.binary)
      
      for (i in 1:nresp){
        for (j in 1:nvotes.binary){
          if (onedim.classification.errors[i,j] == 1) correct.class.binary[i,j] <- 0
          if (onedim.classification.errors[i,j] == 0) correct.class.binary[i,j] <- 1
          if (votemat.binary[i,j] == 9) correct.class.binary[i,j] <- NA
        }}
      
      correct.class.scale <- matrix(NA,nrow=nresp,ncol=nvotes)
      
      for (i in 1:nresp){
        if (is.na(correct.class.binary[i,1])) correct.class.scale[i,1] <- NA
        if (!is.na(correct.class.binary[i,1]) & all(correct.class.binary[i,1:(ncat[1]-1)]==TRUE)) correct.class.scale[i,1] <- 1
        if (!is.na(correct.class.binary[i,1]) & any(correct.class.binary[i,1:(ncat[1]-1)]==FALSE)) correct.class.scale[i,1] <- 0
      }
      
      for (i in 1:nresp){
        for (j in 2:nvotes){
          if (is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)])) correct.class.scale[i,j] <- NA
          if (!is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)]) &
              all(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]==TRUE)) correct.class.scale[i,j] <- 1
          if (!is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)]) &
              any(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]==FALSE)) correct.class.scale[i,j] <- 0
        }}
      
      correct.scale.respondents <- matrix(NA,nrow=nresp,ncol=2)
      
      colnames(correct.scale.respondents) <- c("wrongScale","correctScale")
      
      for (i in 1:nresp){
        correct.scale.respondents[i,] <- table(correct.class.scale[i,])
      }
      
      correct.scale.issues <- matrix(NA,nrow=nvotes,ncol=6)
      
      colnames(correct.scale.issues) <- c("wrongScale","correctScale","modalChoice","nChoices","errorsNull","PREScale")
      
      for (j in 1:nvotes){
        correct.scale.issues[j,1:2] <- table(correct.class.scale[,j])
        correct.scale.issues[j,3] <- max(table(votemat[,j]))
      }
      
      correct.scale.issues[,4] <- correct.scale.issues[,1] + correct.scale.issues[,2]
      correct.scale.issues[,5] <- correct.scale.issues[,4] - correct.scale.issues[,3]
      correct.scale.issues[,6] <- (correct.scale.issues[,5] - correct.scale.issues[,1]) / correct.scale.issues[,5]
      
      issues.unique <- cbind(issues.unique,correct.scale.issues)
      colnames(issues.unique) <- c("correctPerc","PRE","wrongScale","correctScale","modalChoice","nChoices","errorsNull","PREScale")
      
      tempResp <- cbind(correct.scale.respondents,(correct.scale.respondents[,2]/(correct.scale.respondents[,1] + correct.scale.respondents[,2])))
      
      respondents.2 <- matrix(NA, nrow=hr.binary$n, ncol=3)
      respondents.2[use.respondents,] <- tempResp
      
      respondents <- cbind(respondents,respondents.2)
      colnames(respondents) <- c("rank","correctYea","wrongYea","wrongNay","correctNay","volume","coord1D","wrongScale","correctScale","percent.correctScale")
      
      overall.correctClassScale <- sum(issues.unique[,4]) / sum(issues.unique[,6])
      overall.APREscale <- (sum(issues.unique[,7]) - sum(issues.unique[,3])) / (sum(issues.unique[,7]))
      
      fits <- c(fits,overall.correctClassScale,overall.APREscale)
      
      oocObject <- list(respondents = respondents,
                        issues = issues,
                        issues.unique = issues.unique,
                        fits = fits,
                        oc.result.binary = result.binary,
                        votemat.binary = votemat.binary,
                        errorcount = NULL)
      
    }
    
    ##################################################################  
    ## Running Ordered Optimal Classification (Multiple Dimensions) ##
    ##################################################################
    if (ndim>=2){
      
      cat("\n\tRunning Ordered Optimal Classification...\n\n")
      start <- proc.time()
      
      # XMAT = starting values for the respondent coordinates
      
      XMAT <- result.binary$legislators[use.respondents,grep("coord", colnames(result.binary$legislators))]
      XMAT[is.na(XMAT)] <- runif(sum(is.na(XMAT)), -1, 1)
      
      # ZVEC = starting values for the normal vectors
      
      # Random Starts for the Normal Vectors
      Nj <- hitandrun::hypersphere.sample(ndim, nvotes)
      for (j in 1:nvotes){
        if (Nj[j,1] < 0) Nj[j,] <- -1 * Nj[j,]
      }
      
      ZVEC <- matrix(rep(Nj[1,], ncat[1]-1), ncol=ndim, byrow=TRUE)
      for (j in 2:nvotes){
        ZVEC <- rbind(ZVEC, matrix(rep(Nj[j,], (ncat[j]-1)), ncol=ndim, byrow=TRUE))
      }
      
      # LDATA = roll call data
      
      LDATA <- votemat.binary
      nummembers <- nresp
      np <- nresp
      npzz <- np
      numvotes <- ncol(votemat.binary)
      nrcall <- ncol(votemat.binary)
      nv <- nrcall # nrcall = SUM_j=1,q{ncat_j-1}
      ns <- ndim
      ndual <- ndim * (nummembers + numvotes) + 111
      # ndual <- 2 * (nummembers + numvotes) + 111
      LLEGERR <- rep(0,np*ns)
      dim(LLEGERR) <- c(np,ns)
      LERROR <- rep(0,np*nrcall)
      dim(LERROR) <- c(np,nrcall)
      MCUTS <- rep(0,nrcall*2)
      dim(MCUTS) <- c(nrcall,2)
      
      X3 <- rep(0,np*ns)
      dim(X3) <- c(np,ns)
      
      WS <- rep(0,ndual)
      dim(WS) <- c(ndual,1)
      LLV <- rep(0,ndual)
      dim(LLV) <- c(ndual,1)
      LLVB <- rep(0,ndual)
      dim(LLVB) <- c(ndual,1)
      LLE <- rep(0,ndual)
      dim(LLE) <- c(ndual,1)
      LLEB <- rep(0,ndual)
      dim(LLEB) <- c(ndual,1)
      ZS <- rep(0,ndual)
      dim(ZS) <- c(ndual,1)
      XJCH <- rep(0,25)
      dim(XJCH) <- c(25,1)
      XJEH <- rep(0,25)
      dim(XJEH) <- c(25,1)
      XJCL <- rep(0,25)
      dim(XJCL) <- c(25,1)
      XJEL <- rep(0,25)
      dim(XJEL) <- c(25,1)
      errorcount <- rep(0,iter*6)
      dim(errorcount) <- c(iter,6)
      errorcounts <- rep(0,iter*6)
      dim(errorcounts) <- c(iter,6)
      
      YSS <- NULL
      MM <- NULL
      XXX <- NULL
      MVOTE <- NULL
      WSSAVE <- NULL
      
      #  *****************************************************************
      #
      #  GLOBAL LOOP
      #
      #  *****************************************************************
      
      for (iiii in 1:iter){
        
        #  BEGIN JANICE LOOP
        kerrors <- 0
        kchoices <- 0
        kcut <- 1
        lcut <- 6
        for (jx in 1:nrcall){
          for (i in 1:np){
            
            sum=0.0	
            for (k in 1:ns){
              sum <- sum+XMAT[i,k]*ZVEC[jx,k]
            }
            
            #  SAVE PROJECTION VECTORS -- RESPONDENT BY ROLL CALL MATRIX
            
            YSS[i]=sum
            XXX[i]=sum
            MM[i]=LDATA[i,jx]
            if(LDATA[i,jx]==0)MM[i]=9
          }
          
          #  END OF RESPONDENT PROJECTION LOOP
          
          #  SORT PROJECTION VECTOR (Y-HAT)
          
          LLL <- order(XXX)
          YSS <- sort(XXX)
          MVOTE <- MM[LLL]
          
          
          #  FIRST FIND MULTIPLE CUTTING POINTS CORRESPONDING TO EACH NORMAL VECTORS -- THIS STEMS FROM THE STACKING
          #  BY ISSUE EXPLAINED ABOVE
          
          ivot <- jx
          jch <- 0
          jeh <- 0
          jcl <- 0
          jel <- 0
          irotc <- 0
          rescutpoint <- find_cutpoint(nummembers,numvotes,
                                       npzz,nv,np,nrcall,ns,ndual,ivot,XMAT,YSS,MVOTE,WS,
                                       LLV,LLVB,LLE,LLEB,
                                       LERROR,ZS,jch,jeh,jcl,jel,irotc,kcut,lcut,
                                       LLL,XJCH,XJEH,XJCL,XJEL)
          
          #  STORE DIRECTIONALITY OF ROLL CALL
          
          WSSAVE[jx] <- rescutpoint[[13]][jx]
          jch <- rescutpoint[[20]]
          jeh <- rescutpoint[[21]]
          jcl <- rescutpoint[[22]]
          jel <- rescutpoint[[23]]
          kerrors <- kerrors+jeh+jel
          kchoices <- kchoices+jch+jeh+jcl+jel
          kcut <- rescutpoint[[25]]
          lcut <- rescutpoint[[26]]
          MCUTS[jx,1] <- kcut
          MCUTS[jx,2] <- lcut
        }
        
        # END OF ROLL CALL LOOP FOR JANICE
        
        errorcount[iiii,1] <- kerrors
        WS <- WSSAVE
        jjj <- 1
        ltotal <- 0
        mwrong <- 0
        iprint <- 0
        respoly <- find_polytope(nummembers,numvotes,
                                 jjj,np,nrcall,ns,ndual,XMAT,LLEGERR,
                                 ZVEC,WS,MCUTS,LERROR,ltotal,mwrong,
                                 LDATA,iprint)
        
        #  DATA STORED FORTRAN STYLE -- STACKED BY COLUMNS
        
        errorcount[iiii,2] <- respoly[[15]]
        dim(respoly[[8]]) <- c(np,ns)
        X3 <- respoly[[8]]
        r1 <- cor(X3[,1],XMAT[,1])
        r2 <- cor(X3[,2],XMAT[,2])
        errorcount[iiii,3] <- r1
        errorcount[iiii,4] <- r2
        XMAT <- X3
        
        cat("\t\tGetting respondent coordinates...\n")
        
        #  ROLL CALL LOOP (CALCULATE NORMAL VECTORS VIA oprobit,
        #  svm.reg, svm.class, OR krls)
        
        XMATTEMP <- XMAT
        ZVEC_regress <- rep(0,ndim*nvotes)
        dim(ZVEC_regress) <- c(nvotes,ndim)
        
        for (j in 1:nvotes){
          
          # Estimate and Normal Vector
          nvest <- compute.nv(xmat=XMATTEMP, resp=votemat[,j], nv.method=nv.method, ...)
          
          # Store the Output
          for (q in 1:ndim) ZVEC_regress[j,q] <- nvest[q]
          
        }
        
        ZVECpost <- matrix(rep(ZVEC_regress[1,], ncat[1]-1), ncol=ndim, byrow=TRUE)
        for (j in 2:nvotes){
          ZVECpost <- rbind(ZVECpost, matrix(rep(ZVEC_regress[j,], (ncat[j]-1)), ncol=ndim, byrow=TRUE))
        }
        for (j in 1:nrcall){
          if (ZVECpost[j,1] < 0) ZVECpost[j,] <- -1 * ZVECpost[j,]
        }
        r1 <- cor(ZVEC[,1],ZVECpost[,1])
        r2 <- cor(ZVEC[,2],ZVECpost[,2])
        errorcount[iiii,5] <- r1
        errorcount[iiii,6] <- r2
        
        ZVEC <- matrix(rep(ZVEC_regress[1,], ncat[1]-1), ncol=ndim, byrow=TRUE)
        for (j in 2:nvotes){
          ZVEC <- rbind(ZVEC, matrix(rep(ZVEC_regress[j,], (ncat[j]-1)), ncol=ndim, byrow=TRUE))
        }
        
        for (j in 1:nrcall){
          if (ZVEC[j,1] < 0) ZVEC[j,] <- -1 * ZVEC[j,]
        }
        
        cat("\t\tCalculating normal vectors...\n")
        
      }
      #  END OF GLOBAL LOOP
      
      #  RUN POLYTOPE SEARCH ONE FINAL TIME
      respoly <- find_polytope(nummembers,numvotes,
                               jjj,np,nrcall,ns,ndual,XMAT,LLEGERR,
                               ZVEC,WS,MCUTS,LERROR,ltotal,mwrong,
                               LDATA,iprint)
      
      dim(respoly[[8]]) <- c(np,ns)
      X3 <- respoly[[8]]
      XMAT <- X3
      totalerrors <- respoly[[15]]
      
      polarity <- matrix(NA,nrow=nresp,ncol=sum(ncat)-nvotes)
      correctYea <- NA
      wrongYea <- NA
      wrongNay <- NA
      correctNay <- NA
      kchoices <- NA
      kminority <- NA
      voteerrors <- NA
      PRE <- NA
      respondent.class.correctYea <- matrix(0, nrow=nresp, ncol=(sum(ncat)-nvotes))
      respondent.class.wrongYea <- matrix(0, nrow=nresp, ncol=(sum(ncat)-nvotes))
      respondent.class.wrongNay <- matrix(0, nrow=nresp, ncol=(sum(ncat)-nvotes))
      respondent.class.correctNay <- matrix(0, nrow=nresp, ncol=(sum(ncat)-nvotes))
      
      for (j in 1:(sum(ncat)-nvotes)){
        ivote <- votemat.binary[,j]
        #polarity <- XMAT[,1]*ZVEC[j,1] + XMAT[,2]*ZVEC[j,2] - WS[j]
        polarity <- as.vector(XMAT%*%ZVEC[j,]) - WS[j]
        errors1 <- ivote==1 & polarity >= 0
        errors2 <- ivote==6 & polarity <= 0
        errors3 <- ivote==1 & polarity <= 0
        errors4 <- ivote==6 & polarity >= 0
        kerrors1 <- ifelse(is.na(errors1),9,errors1)
        kerrors2 <- ifelse(is.na(errors2),9,errors2)
        kerrors3 <- ifelse(is.na(errors3),9,errors3)
        kerrors4 <- ifelse(is.na(errors4),9,errors4)
        kerrors12 <- kerrors34 <- 0
        kerrors12 <- sum(kerrors1==1)+sum(kerrors2==1)
        kerrors34 <- sum(kerrors3==1)+sum(kerrors4==1)
        if (kerrors12 >= kerrors34){
          yeaerror <- errors3
          nayerror <- errors4
        }
        if (kerrors12 < kerrors34){
          yeaerror <- errors1
          nayerror <- errors2
        }
        kerrorsmin <- min(kerrors12,kerrors34)
        
        respondent.class.correctYea[,j] <- as.numeric(ivote==1 & !yeaerror & !nayerror)
        respondent.class.wrongYea[,j] <- as.numeric(ivote==6 & nayerror)
        respondent.class.wrongNay[,j] <- as.numeric(ivote==1 & yeaerror)
        respondent.class.correctNay[,j] <- as.numeric(ivote==6 & !yeaerror & !nayerror)
        
        kpyea <- sum(ivote==1)
        kpnay <- sum(ivote==6)
        kpyeaerror <- sum(yeaerror)
        kpnayerror <- sum(nayerror)
        correctYea[j] <- kpyea - kpyeaerror
        wrongYea[j] <- kpnayerror
        wrongNay[j] <- kpyeaerror
        correctNay[j] <- kpnay - kpnayerror
        kchoices[j] <- kpyea+kpnay
        kminority[j] <- min(kpyea,kpnay)
        voteerrors[j] <- kerrorsmin
        PRE[j] <- (min(kpyea,kpnay) - kerrorsmin)/(min(kpyea,kpnay))
      }
      
      totalerrors2 <- sum(voteerrors)
      totalcorrect2 <- sum(correctYea) + sum(correctNay)
      APRE <- ((sum(kminority) - sum(voteerrors)) / sum(kminority))
      totalCorrectClass <- totalcorrect2 / (totalcorrect2 + totalerrors2)
      
      normVectorkD <- ZVEC
      midpoints <- WS
      
      issues <- cbind(correctYea,wrongYea,wrongNay,correctNay,PRE,normVectorkD,midpoints)
      colnames(issues) <- c("correctYea","wrongYea","wrongNay","correctNay","PRE",paste("normVector",1:ndim,"D",sep=""),"midpoints")
      
      respondent.correctYea <- rowSums(respondent.class.correctYea)
      respondent.wrongYea <- rowSums(respondent.class.wrongYea)
      respondent.wrongNay <- rowSums(respondent.class.wrongNay)
      respondent.correctNay <- rowSums(respondent.class.correctNay)
      
      correct.class.binary <- matrix(NA,nrow=nresp,ncol=ncol(respondent.class.correctYea))
      
      for (i in 1:nresp){
        for (j in 1:ncol(respondent.class.correctYea)){
          
          if (respondent.class.correctYea[i,j] == 1) correct.class.binary[i,j] <- 1
          if (respondent.class.wrongYea[i,j] == 1) correct.class.binary[i,j] <- 0
          if (respondent.class.wrongNay[i,j] == 1) correct.class.binary[i,j] <- 0
          if (respondent.class.correctNay[i,j] == 1) correct.class.binary[i,j] <- 1
          
        }}
      
      correct.class.scale <- matrix(NA,nrow=nresp,ncol=nvotes)
      
      for (i in 1:nresp){
        if (is.na(correct.class.binary[i,1])) correct.class.scale[i,1] <- NA
        if (!is.na(correct.class.binary[i,1]) & all(correct.class.binary[i,1:(ncat[1]-1)]==TRUE)) correct.class.scale[i,1] <- 1
        if (!is.na(correct.class.binary[i,1]) & any(correct.class.binary[i,1:(ncat[1]-1)]==FALSE)) correct.class.scale[i,1] <- 0
      }
      
      for (i in 1:nresp){
        for (j in 2:nvotes){
          if (is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)])) correct.class.scale[i,j] <- NA
          if (!is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)]) &
              all(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]==TRUE)) correct.class.scale[i,j] <- 1
          if (!is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)]) &
              any(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]==FALSE)) correct.class.scale[i,j] <- 0
        }}
      
      correct.scale.respondents <- matrix(NA,nrow=nresp,ncol=2)
      
      colnames(correct.scale.respondents) <- c("wrongScale","correctScale")
      
      for (i in 1:nresp){
        correct.scale.respondents[i,] <- table(correct.class.scale[i,])
      }
      
      correct.scale.issues <- matrix(NA,nrow=nvotes,ncol=6)
      
      colnames(correct.scale.issues) <- c("wrongScale","correctScale","modalChoice","nChoices","errorsNull","PREScale")
      
      for (j in 1:nvotes){
        correct.scale.issues[j,1:2] <- table(correct.class.scale[,j])
        correct.scale.issues[j,3] <- max(table(votemat[,j]))
      }
      
      correct.scale.issues[,4] <- correct.scale.issues[,1] + correct.scale.issues[,2]
      correct.scale.issues[,5] <- correct.scale.issues[,4] - correct.scale.issues[,3]
      correct.scale.issues[,6] <- (correct.scale.issues[,5] - correct.scale.issues[,1]) / correct.scale.issues[,5]
      
      tempResp <- cbind(respondent.correctYea,respondent.wrongYea,respondent.wrongNay,respondent.correctNay,XMAT,
                        correct.scale.respondents,(correct.scale.respondents[,2]/(correct.scale.respondents[,1] + correct.scale.respondents[,2])))
      
      respondents <- matrix(NA, nrow=hr.binary$n, ncol=(7+ndim))
      respondents[use.respondents,] <- tempResp
      
      colnames(respondents) <- c("correctYea","wrongYea","wrongNay","correctNay",paste("coord",1:ndim,"D",sep=""),"wrongScale","correctScale","percent.correctScale")
      
      # COMPRESS ROLL CALL MATRIX
      
      unique.normVectorkD <- unique(normVectorkD)
      unique.NV1 <- unique.normVectorkD[,1]
      unique.NV2 <- unique.normVectorkD[,2]
      minorityvote <- apply(cbind(issues[,1] + issues[,3], issues[,2] + issues[,4]),1,min)
      classerrors <- issues[,2] + issues[,3]
      correctclass <- issues[,1] + issues[,4]
      
      uniqueCLASSPERC <- sum(correctclass[1:(ncat[1]-1)]) / sum(correctclass[1:(ncat[1]-1)] + classerrors[1:(ncat[1]-1)])
      for (j in 2:nvotes){
        uniqueCLASSPERC[j] <- sum(correctclass[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]) / sum(correctclass[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])] + classerrors[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])])
      }
      
      uniquePRE <- sum(minorityvote[1:(ncat[1]-1)] - classerrors[1:(ncat[1]-1)]) / sum(minorityvote[1:(ncat[1]-1)])
      for (j in 2:nvotes){
        uniquePRE[j] <- sum(minorityvote[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])] - classerrors[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]) / sum(minorityvote[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])])
      }
      
      b <- c(1.0,0)
      theta <- NULL
      for (j in 1:nvotes){
        a <- c(unique.NV1[j], unique.NV2[j])
        theta[j] <- aspace::acos_d(sum(a*b) / (sqrt(sum(a*a))*sqrt(sum(b*b))))
        if (a[2] < 0) theta[j] <- -1 * theta[j]
      }
      
      issues.unique <- cbind(uniqueCLASSPERC,uniquePRE,unique.normVectorkD,theta,correct.scale.issues)
      colnames(issues.unique) <- c("correctPerc","PRE",paste("normVector",1:ndim,"D",sep=""),"normVectorAngle2D",
                                   "wrongScale","correctScale","modalChoice","nChoices","errorsNull","PREScale")
      
      overall.correctClassScale <- sum(issues.unique[,"correctScale"]) / sum(issues.unique[,"nChoices"])
      overall.APREscale <- (sum(issues.unique[,"errorsNull"]) - sum(issues.unique[,"wrongScale"])) / (sum(issues.unique[,"errorsNull"]))
      
      fits <- c(totalCorrectClass,APRE,overall.correctClassScale,overall.APREscale)
      
      # <errorcount>
      #   First column = errors after cutpoint search
      #   Second column = errors after polytope search
      #   Third column = Correlation respondent coordinates 1st dim
      #   Fourth column = Correlation respondent coordinates 2nd dim
      #   Fifth column = Correlation normal vectors 1st dim
      #   Sixth column = Correlation normal vectors 2nd dim
      
      oocObject <- list(respondents = respondents,
                        issues = issues,
                        issues.unique = issues.unique,
                        fits = fits,
                        oc.result.binary = result.binary,
                        errorcount = errorcount,
                        votemat.binary = votemat.binary)
      
      cat("\n\nOrdered Optimal Classification estimation completed successfully.")
      cat("\nOrdered Optimal Classification took", (proc.time() - start)[3], "seconds to execute.\n\n")
      
    }
  
  } else if (binary.method=="even") {
    
    # the matrix with binary choices recoded to be as evenly-balanced as possible
    votemat.binary <- apply(votemat, 2, binary.balanced)
    
    hr.binary <- rollcall(votemat.binary, yea = 1, 
                               nay = 6, missing = 9, notInLegis = NA)
    result.binary <- oc(hr.binary, dims = dims, minvotes = minvotes, 
                             lop = lop, polarity = polarity)
    
    
    oocObject <- list(respondents = result.binary$legislators,
                      issues = result.binary$rollcalls,
                      issues.unique = NULL,
                      fits = result.binary$fits,
                      oc.result.binary = result.binary,
                      votemat.binary = votemat.binary,
                      errorcount = NULL)
    
    
  } else {
    
    stop("Invalid 'binary.method' value!")
    
  }
  
  ## Additional Information
  oocObject$dimensions <- ndim
  
  ## Export Output
  class(oocObject) <- "oocflexObject"
  oocObject

}

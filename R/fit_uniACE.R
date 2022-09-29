#' fit_uniACE
#' @description Use OpenMx to quickly fit a univariate ACE model
#' @import OpenMx
#' @param data_1 A n by 2 \code{data.frame} consisting of the group1 kin pairs
#' @param data_2 A n by 2 \code{data.frame} consisting of the group2 kin pairs
#' @param GroupRel A numeric vector specifying two genetic relatedness values of two groups of kin pairs
#' @param GroupR_c A numeric vector specifying two common environment correlation coefficients of two groups of kin pairs
#' @param lbound A logical value indicating if a lower boundary of .0001 will be imposed to the estimated A, C and E components
#' @return Returns a \code{list} with the following:
#' \item{df_nested}{A \code{data.frame} displaying the nested comparison model between ACE, AE, CE, E models}
#' \item{fitACE}{A \code{list} of all model fit information generated from OpenMx}
#' @export

fit_uniACE <- function(data_1, data_2, GroupRel = c(1,.5), GroupR_c = c(1,1), lbound = FALSE){
     # Load Libraries & Options
     #require(OpenMx)
     #require(psych)
     #require(polycor)
     # source("miFunctions.R")
     # # Create Output
     # filename <- "oneACEc"
     # sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

     # ----------------------------------------------------------------------------------------------------------------------
     # PREPARE DATA
     # Load Data
     mzData    <- data_1
     dzData    <- data_2
     R1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1,free = FALSE, values = GroupRel[1], name = "R1")
     R2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1,free = FALSE, values = GroupRel[2], name = "R2")
     r_c1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1,free = FALSE, values = GroupR_c[1], name = "r_c1")
     r_c2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1,free = FALSE, values = GroupR_c[2], name = "r_c2")
     #coeAM <- coe_am

     # covMZ <- cov(mzData, use = "pairwise")
     # covDZ <- cov(dzData, use = "pairwise")
     # #
     # mean(rbind(mzData,dzData)[,1], na.rm = TRUE)

     nv <- 1
     ntv <- 2
     selVars1   <- colnames(mzData)
     selVars2   <- colnames(dzData)

     #start values
     svBe <- .01
     svMu <- 0
     svVa <- .2
     svVe <- .5
     V <- NULL
     VA <- NULL
     VC <- NULL
     VE <- NULL
     cDZ <- NULL
     cMZ <- NULL
     #variance matrix

     if(lbound == TRUE){
          covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVa, lbound = .0001, labels = "VA11", name = "VA")
          covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVa, lbound = .0001, labels = "VC11", name = "VC")
          covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVe, lbound = .0001, labels = "VE11", name = "VE")
     }else{
          covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVa, labels = "VA11", name = "VA")
          covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVa, labels = "VC11", name = "VC")
          covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVe, labels = "VE11", name = "VE")
     }


     #expected variance matrix
     covP <- mxAlgebra(expression = VA+VC+VE, name = "V")
     covMZ <- mxAlgebra(expression = R1*VA+r_c1*VC, name = "cMZ")
     covDZ <- mxAlgebra(expression = R2%x%VA+r_c2*VC, name = "cDZ")
     expCovMz <- mxAlgebra(expression = rbind(cbind(V,cMZ), cbind(t(cMZ),V)), name = "expCovMz")
     expCovDz <- mxAlgebra(expression = rbind(cbind(V,cDZ), cbind(t(cDZ),V)), name = "expCovDz")

     #create data
     dataMZ       <- mxData( observed=mzData, type="raw" )
     dataDZ       <- mxData( observed=dzData, type="raw" )

     # Mean Matrix
     intercept <- mxMatrix(type = "Full", nrow= 1 , ncol = ntv, free = TRUE, values = 0, labels = "interC", name = "intercept")
     expMean <- mxAlgebra(expression = 1*intercept , name = "expMean")

     # Create expectation objects
     expMZ <- mxExpectationNormal(covariance = "expCovMz", means ="expMean", dimnames = selVars1)
     expDZ <- mxExpectationNormal(covariance = "expCovDz", means ="expMean", dimnames = selVars2)
     funML <- mxFitFunctionML()

     #Create models
     pars <- list(intercept, covA, covC, covE, covP)
     modelMZ <- mxModel(pars, expMean,covMZ,expCovMz,dataMZ,expMZ,funML,R1,r_c1,name = "MZ")
     #MZfit <- mxRun(modelMZ, intervals = TRUE)
     #summary(MZfit)
     modelDZ <- mxModel(pars, expMean,covDZ,expCovDz,dataDZ,expDZ,funML,R2,r_c2,name = "DZ")
     #DZfit <- mxRun(modelDZ, intervals = TRUE)
     #summary(DZfit)

     multi <- mxFitFunctionMultigroup(c("MZ","DZ"))


     #Algebra for Variance components
     rowUS <- rep("US",nv)
     colUS <- rep(c("VA","VC","VE","SA","SC","SE"),each = nv)
     estUS <- mxAlgebra(expression = cbind(VA,VC,VE,VA/V,VC/V,VE/V), name = "US", dimnames = list(rowUS,colUS))

     #CI
     ciACE <- mxCI("US[1,1:6]")
     modelACE <- mxModel("oneACEvc_1cov", pars, modelMZ, modelDZ, multi,estUS,ciACE)
     fitACE <- mxRun(modelACE, intervals = TRUE, silent = TRUE)
     sumACE <-summary(fitACE)
     #sumACE

     # ----------------------------------------------------------------------------------------------------------------------
     # RUN SUBMODELS
     # Run AE model
     modelAE <- mxModel( modelACE, name="oneAEvc" )
     modelAE <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
     fitAE <- mxRun( modelAE, intervals=T, silent = TRUE )
     #fitGofs(fitAE); fitEstCis(fitAE)
     # Run CE model
     modelCE <- mxModel( modelACE, name="oneCEvc" )
     modelCE <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
     modelCE <- omxSetParameters( modelCE, labels=c("VE11","VC11"), free=TRUE, values=.6 )
     fitCE <- mxRun( modelCE, intervals=T, silent = TRUE )
     #fitGofs(fitCE); fitEstCis(fitCE)
     # Run E model
     modelE <- mxModel( modelACE, name="oneEvc" )
     modelE <- omxSetParameters( modelE, labels=c("VA11","VC11"), free=FALSE, values=0 )
     fitE <- mxRun( modelE, intervals=T, silent = TRUE )
     #fitGofs(fitE); fitEstCis(fitE)
     # Print Comparative Fit Statistics
     df_nested <- mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE) )
     #(rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result,fitE$US$result),4)
     l.modeloutput <- list(nest = df_nested,summary = sumACE)
     return(l.modeloutput)
}

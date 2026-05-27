#' fit_OrdACE
#' @description Use OpenMx to quickly fit a univariate Ordinal ACE model
#' @importFrom OpenMx mxMatrix mxAlgebra mxData mxExpectationNormal mxFitFunctionML mxModel mxRun mxCI mxFitFunctionMultigroup omxSetParameters mxCompare
#' @importFrom polycor hetcor
#' @param data_1 A n by 2 \code{data.frame} consisting of the group1 kin pairs
#' @param data_2 A n by 2 \code{data.frame} consisting of the group2 kin pairs
#' @param GroupRel A numeric vector specifying two genetic relatedness values of two groups of kin pairs
#' @param GroupR_c A numeric vector specifying two common environment correlation coefficients of two groups of kin pairs
#' @param nth A numerical value specifiying the number of thresholds
#' @param lbound A logical value indicating if a lower boundary of .0001 will be imposed to the estimated A, C and E components
#' @return Returns a \code{list} with the following:
#' \item{df_nested}{A \code{data.frame} displaying the nested comparison model between ACE, AE, CE, E models}
#' \item{fitACE}{A \code{list} of all model fit information generated from OpenMx}
#' @export

fit_OrdACE <- function(data_1, data_2, GroupRel = c(1, .5), GroupR_c = c(1, 1),
                       nth = 4
                       , lbound = FALSE) {
  # Load Libraries & Options
  # require(OpenMx)
  # require(psych)
  # require(polycor)
  # source("miFunctions.R")
  # # Create Output
  # filename <- "oneACEc"
  # sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

  # ----------------------------------------------------------------------------------------------------------------------
  # PREPARE DATA

 # nth <- nth

  # Load Data
  FSData <- data_1 #in the Sim_Fit2.R function it already assigns the groups and the variable names
  HSData <- data_2
  FSDataF <-  mxFactor( x=FSData, levels=c(0:nth) )
  HSDataF <- mxFactor( x=HSData, levels=c(0:nth) )

  vars <- 'Ord_' # I don't know what to do with this just yet - this is the list of variables
  nv <- 1 # number of variables
  ntv <- nv*2 # number of total variables
  selVars <- c("Ord_1", "Ord_2") #paste(vars,c(rep(1,nv),rep(2,nv)),sep="")


  R1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, values = GroupRel[1], name = "R1")
  R2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, values = GroupRel[2], name = "R2")
  r_c1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, values = GroupR_c[1], name = "r_c1")
  r_c2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, values = GroupR_c[2], name = "r_c2")


  #Descriptives
  sapply(FSData,table)
  sapply(HSData,table)
  hetcor(FSData)$cor
  hetcor(HSData)$cor

  # coeAM <- coe_am

  # covMZ <- cov(mzData, use = "pairwise")
  # covDZ <- cov(dzData, use = "pairwise")
  # #
  # mean(rbind(mzData,dzData)[,1], na.rm = TRUE)

  nv <- 1
  ntv <- 2
  selVars1 <- colnames(FSData)
  selVars2 <- colnames(HSData)

  # start values
  svLTh <- 0.01 # start value for first threshold
  svITh <- 1 # start value for increments
  svTh <- matrix(rep(c(svLTh,(rep(svITh,nth-1)))),nrow=nth,ncol=nv) # start value for thresholds
  lbTh <- matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv) # lower bounds for thresholds
  svPa <- .2 # start value for path coefficient
  svPc <- .3
  svPe <- .4 #start value for the path coefficient e

  # variance matrix

  # if (lbound == TRUE) {
  #   covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, values = svVa, lbound = .0001, labels = "VA11", name = "VA")
  #   covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, values = svVa, lbound = .0001, labels = "VC11", name = "VC")
  #   covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, values = svVe, lbound = .0001, labels = "VE11", name = "VE")
  # } else {
  #   covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, values = svVa, labels = "VA11", name = "VA")
  #   covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, values = svVa, labels = "VC11", name = "VC")
  #   covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, values = svVe, labels = "VE11", name = "VE")
  # }

  #PREPARE MODEL
  # Create Algebra for expected Mean & Threshold Matrices

  meanG <- mxMatrix( type="Zero", nrow=1, ncol=ntv, name="meanG" )

  thinG <- mxMatrix( type="Full", nrow=nth, ncol=ntv, free=TRUE, values=svTh, lbound=lbTh, labels=labTh("th",vars,nth), name="thinG")

  inc <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=FALSE, values=1, name="inc" )

  threG <- mxAlgebra( expression= inc %*% thinG, name="threG" )

  #Create matrices for variance components
  covA <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11", name="VA" )
  covC <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VC11", name="VC" )
  covE <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VE11", name="VE") #, lbound = 0.0001 ) makes it so that E isn't negative

  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covP <- mxAlgebra( expression= VA+VC+VE, name="V" )
  covFS <- mxAlgebra( expression= 0.5%x%VA+VC, name="cFS" )
  covHS <- mxAlgebra( expression= 0.25%x%VA+ VC, name="cHS" )
  expCovFS <- mxAlgebra( expression= rbind( cbind(V, cFS), cbind(t(cFS), V)), name="expCovFS" )
  expCovHS <- mxAlgebra( expression= rbind( cbind(V, cHS), cbind(t(cHS), V)), name="expCovHS" )

  # Constrain Variance of Binary Variables
  var1 <- mxConstraint( expression=diag2vec(V)==1, name="Var1" )

  # Create Data Objects for Multiple Groups
  dataFS <- mxData( observed=FSDataF, type="raw" )
  dataHS <- mxData( observed=HSDataF, type="raw" )

  # Create Expectation Objects for Multiple Groups
  expFS <- mxExpectationNormal( covariance="expCovFS", means="meanG", dimnames=selVars, thresholds="threG" )
  expHS <- mxExpectationNormal( covariance="expCovHS", means="meanG", dimnames=selVars, thresholds="threG" )
  funML <- mxFitFunctionML()

  # Create Model Objects for Multiple Groups
  pars <- list(meanG, thinG,inc, threG,covA, covC, covE, covP)
  modelFS <- mxModel(pars, covFS, expCovFS, dataFS, expFS, funML, name ="FS" )
  modelHS <- mxModel(pars, covHS, expCovHS, dataHS, expHS, funML, name ="HS" )
  multi <- mxFitFunctionMultigroup( c("FS","HS") )

  # Create Algebra for Unstandardized and Standardized Variance Components
  rowUS <- rep('US',nv)
  colUS <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
  estUS <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )

  # Create Confidence Interval Objects
  ciACE <- mxCI( "US[1,1:3]" )

  # Build Model with Confidence Intervals
  modelACE <- mxModel( "oneACEvo", pars, var1, modelFS, modelHS, multi, estUS, ciACE)

  #---------------------------------------------------------------------------------------
  # RUN MODEL
  slsqp <- mxOption(NULL,"Default optimizer","SLSQP")
  # Run ACE Model
  fitACE <- mxTryHardOrdinal( modelACE, intervals=TRUE )
  sumACE <- summary( fitACE )

  # Compare with Saturated Model

  #if saturated model fitted in same session
  #mxCompare( fitSAT, fitACE )
  #if saturated model prior to genetic model
  #lrtSAT(fitACE,4207.7738,1762)

  # Print Goodness-of-fit Statistics & Parameter Estimates
  fitGofs(fitACE)
  fitEstCis(fitACE)

  # ----------------------------------------------------------------------------------------------------------------------
  # RUN SUBMODELS
  # Run AE model
  modelAE <- mxModel( fitACE, name="oneAEvo" )
  modelAE <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
  fitAE <- mxRun( modelAE, intervals=T )
 fitGofs(fitAE); fitEstCis(fitAE)

  # Run CE model
  modelCE <- mxModel( fitACE, name="oneCEvo" )
  modelCE <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
  modelCE <- omxSetParameters( modelCE, labels=c("VE11","VC11"), free=TRUE, values=.6 )
  fitCE <- mxRun( modelCE, intervals=TRUE )
 fitGofs(fitCE); fitEstCis(fitCE)

  # Run E model
 # modelE <- mxModel( fitAE, name="oneEvo" )
 # modelE <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
 # fitE <- mxRun( modelE, intervals=T )
 # fitGofs(fitE); fitEstCis(fitE)

  # Print Comparative Fit Statistics
  mxCompare( fitACE, nested <- list(fitAE, fitCE) ) #  fitE commented out for now
  round(rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result),4) #,fitE$US$result

  #in response to one of my warnings about the optimizer having a status code 5 for modelETO
  #mxCheckIdentification(modelETO, details=TRUE, nrows=2, exhaustive=FALSE, silent=FALSE)

  sumACE

  #this gives me a table with two rows and three columns - row 1 is unstandardized ACE, row 2 is standardized.
  #I got this from https://openmx.ssri.psu.edu/docs/OpenMx/2.5.1/GeneticEpi_Path.html

  # Generate & Print Output
  # additive genetic variance, a^2
  A  <- mxEval(VA11, fitACE)
  # shared environmental variance, c^2
  C  <- mxEval(VC11, fitACE)
  # unique environmental variance, e^2
  E  <- mxEval(VE11, fitACE)
  # total variance
  V  <- (A+C+E)
  # standardized A
  a2 <- A/V
  # standardized C
  c2 <- C/V
  # standardized E
  e2 <- E/V
  # table of estimates
  estACE <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
  # likelihood of ACE model
  LL_ACE <- mxEval(fitfunction, fitACE)


   # Print Comparative Fit Statistics

   df_nested <- mxCompare(fitACE, nested <- list(fitAE, fitCE ))#fitE) commented out for now
  # (rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result,fitE$US$result),4)
  l.modeloutput <- list(nest = df_nested, summary = sumACE)
  return(l.modeloutput)
}



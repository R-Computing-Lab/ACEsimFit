#' Sim_Fit
#' @description A function to simulate a set of kin pair data and fit them with ACE models. Can be helpful with checking model performance for a given parameter setting.
#' @param GroupNames A character vector specifying two names of the simulated kin pairs
#' @param GroupSizes A numeric vector specifying two group sizes indicating the amount of kin pairs in respective group.
#' @param nIter A numeric value specifying the number of iteration you want to run given the parameters assigned (i.e. the number of model fitting results you want to get)
#' @param SSeed An integer specifying the starting seed of the random number. This parameter will make sure the simulated results are replicable across time
#' @param GroupRel A numeric vector specifying two genetic relatedness values of the simulated kin pairs
#' @param GroupR_c A numeric vector specifying two common environment correlation coefficients of the simulated kin pairs
#' @param mu A numeric vector specifying two mean values for the generated variable of the kin pairs
#' @param ace1 A numeric vector specifying three variance components under an ACE (additive genetics, common environment, unique environment) structure for group1
#' @param ace2 A numeric vector specifying three variance components under an ACE (additive genetics, common environment, unique environment) structure for group2
#' @param ifComb A logical value specifying the approach to achieve the required genetic relatedness value. \code{TRUE} = using combination approach. \code{FALSE} = using direct approach. (See function description for a detailed explanation of two approaches.)
#' @param lbound A logical value indicating if a lower boundary of .0001 will be imposed to the estimated A, C and E components
#' @param saveRaw A logical value specifying if the raw simulated data should be saved in the output list
#' @return Returns a two-level \code{list}. Level-one is the number of iterations. Level-two is the model fitting results and raw data (if \code{saveRaw = TRUE}) of the simulated data from the respective iteration. Level-two includes:
#' \item{Results}{A \code{list} including 1) A \code{data.frame} displaying the nested comparison model between ACE, AE, CE, E models and 2) A \code{list} of all model fit information generated from OpenMx}
#' \item{Data}{A \code{data.frame} consists of the simulated raw data}
#' @export

Sim_Fit <- function(
     GroupNames = c("KinPair1","KinPair2"),
     GroupSizes = c(100,100),
     nIter = 100,
     SSeed = 62,
     GroupRel=c(1,.5),
     GroupR_c=c(1,1),
     mu = c(0,0),
     ace1=c(1,1,1),
     ace2=c(1,1,1),
     ifComb=FALSE,
     lbound=FALSE,
     saveRaw=FALSE
){
     l.results <- list()
     for(i in 1: nIter){
          set.seed(SSeed - 1 + i)
          df_temp <- kinsim_double(
               GroupNames = GroupNames,
               GroupSizes = GroupSizes,
               GroupRel=GroupRel,
               GroupR_c=GroupR_c,
               mu = mu,
               ace1=ace1,
               ace2=ace2,
               ifComb=ifComb)
          if(!saveRaw){
               l.results[[i]] <- list(
                    Results =  fit_uniACE(
                         data_1 = df_temp[which(df_temp$GroupName == GroupNames[1]), c("y1","y2")],
                         data_2 = df_temp[which(df_temp$GroupName == GroupNames[2]), c("y1","y2")],
                         GroupRel = GroupRel, GroupR_c = GroupR_c, lbound = lbound),
                    data = NA)
          }else{
               l.results[[i]] <- list(
                    Results =  fit_uniACE(
                         data_1 = df_temp[which(df_temp$GroupName == GroupNames[1]), c("y1","y2")],
                         data_2 = df_temp[which(df_temp$GroupName == GroupNames[2]), c("y1","y2")],
                         GroupRel = GroupRel, GroupR_c = GroupR_c, lbound = lbound),
                    data =  df_temp)
          }
          names(l.results)[[i]] <- paste("Iteration",i, sep = "")
     }
     return(l.results)
}

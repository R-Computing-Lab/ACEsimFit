#' kinsim_double
#' @description The function to generate two groups of univariate kin pair(e.g., both MZ and DZ twins) data using a multivariate norm approach, given the ACE components.
#' @param GroupNames A character vector specifying two names of the simulated kin pairs
#' @param GroupSizes A numeric vector specifying two group sizes indicating the amount of kin pairs in respective group.
#' @param GroupRel A numeric vector specifying two genetic relatedness values of the simulated kin pairs
#' @param GroupR_c A numeric vector specifying two common environment correlation coefficients of the simulated kin pairs
#' @param mu A numeric value specifying the mean for the generated variable
#' @param ace1 A numeric vector specifying three variance components under an ACE (additive genetics, common environment, unique environment) structure for group1
#' @param ace2 A numeric vector specifying three variance components under an ACE (additive genetics, common environment, unique environment) structure for group2
#' @return Returns \code{data.frame} with the following:
#' \item{GroupName}{group name of the kin pairs}
#' \item{R}{level of relatedness for the kin pair}
#' \item{r_c}{level of common envrionment correlation of the kin pairs}
#' \item{id}{id}
#' \item{A1}{Additive genetic component for kin1 of the kin pairs}
#' \item{A2}{Additive genetic component for kin2 of the kin pairs}
#' \item{C1}{shared-environmental component for kin1 of the kin pairs}
#' \item{C2}{shared-environmental component for kin2 of the kin pairs}
#' \item{E1}{non-shared-environmental component for kin1 of the kin pairs}
#' \item{E2}{non-shared-environmental component for kin2 of the kin pairs}
#' \item{y1}{generated variable i for kin1}
#' \item{y2}{generated variable i for kin2}
#' @export

kinsim_double <- function(
     GroupNames = c("KinPair1","KinPair2"),
     GroupSizes = c(100,100),
     GroupRel=c(1,.5),
     GroupR_c=c(1,1),
     mu = 0,
     ace1=c(1,1,1),
     ace2=c(1,1,1)){
     df_N1 <- kinsim_single(
          name = GroupNames[1],
          Rel = GroupRel[1],
          r_c = GroupR_c[1],
          n = GroupSizes[1],
          mu = mu,
          ace = ace1)

     df_N2 <- kinsim_single(
          name = GroupNames[2],
          Rel = GroupRel[2],
          r_c = GroupR_c[2],
          n = GroupSizes[2],
          mu = mu,
          ace = ace2)
     df_final <- rbind(df_N1, df_N2)
     return(df_final)
}

#' kinsim_single
#' @description The function to generate one group of univariate kin pair (e.g., only DZ twins) data using a multivariate norm approach, given the ACE components.
#' @param name Assigned name for the simulated group of kin pairs
#' @param Rel Genetic relatedness of the simulated kin pairs
#' @param r_c Assumed common enviroment correlation
#' @param n The number of generated kin pairs.(n PAIRS of data; The total number of participants is 2n)
#' @param mu The mean for generated variable
#' @param ace Vector of variance components under an ACE (additive genetics, common environment, unique environment) structure
#' @return Returns \code{data.frame} with the following:
#' \item{GroupName}{group name of the kin pairs}
#' \item{R}{level of genetic relatedness for the kin pairs}
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

kinsim_single <- function(
     name = "KinPair1",
     Rel=1,
     r_c = 1,
     n=100,
     mu=0,
     ace=c(1,1,1)){

     sA <- ace[1]^0.5
     sC <- ace[2]^0.5
     sE <- ace[3]^0.5

     S2 <- matrix(c(0,1,
                    1,0),2)
     datalist <- list()


     id <- 1:sum(n)

     A.r <- sA*rmvn(n,
                    sigma = diag(2) + S2*Rel)
     C.r <- sC*rmvn(n,
                    sigma = diag(2) + S2*r_c)
     E.r <- cbind(stats::rnorm(n,
                               sd = sE),
                  stats::rnorm(n,
                               sd = sE))

     y.r <- mu + A.r + C.r + E.r


     r_ <- rep(Rel,n)
     r_c <- rep(r_c,n)

     groupName <- rep(name,n)

     data.r <- data.frame(groupName, r_, r_c, id, A.r,C.r,E.r,y.r)

     names(data.r) <- c("GroupName","R","r_c","id","A1","A2","C1","C2","E1","E2","y1","y2")

     return(data.r)
}

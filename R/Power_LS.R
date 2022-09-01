#' Power_LS
#' @description The function is designed for calculating the power of heritability estimation from ACE models given the parameter settings. Or calculate one of the parameter settings (N,R,h2,c2) given the rest of known parameters.
#' This power calculator is made based on the Least Squares theory and follows the mathematical derivation proposed by Visscher(2004).
#' @import stats
#' @param N1 The number of kin pairs for group1 (amount of PAIRS)
#' @param N2 The number of kin pairs for group2
#' @param power The power of heritability estimation. Specified if you want to return the required sample sizes.
#' @param p_N1 The proportion of kin group1 over the . Required to be specified if the user wants to calculate the N1 and N2 simultaneously.
#' @param h2 The assumed standard heritability value of the target trait. 0 < h2 < 1
#' @param c2 The assumed standard common environmental effects on the target trait. 0 < c2 < 2
#' @param R1 The genetic relatedness of kin pair group1
#' @param R2 The genetic relatedness of kin pair group2
#' @param alpha The type-one error rate for heritability estimation.
#' @return A numeric \code{vector} of power when `N1` and `N2` are both specified.
#' \cr
#' A numeric \code{vector} of `N1` (or `N2`) when `N2` (or `N1`) is specified. A numeric \code{vector} of `N1` and `N2` when `RatioN` is specified.
#' @export


Power_LS <- function(N1,N2,power,p_N1=NULL,h2,c2,R1=1,R2=.5,alpha = .05){
     Za <- qnorm(1-alpha)

     if(missing(power)){
          Zb <- sqrt(h2^2*(abs(R1-R2)^2) / ((1-(R1*h2+c2)^2)^2/N1 + (1-(R2*h2+c2)^2)^2/N2)) - Za
          power_result <- pnorm(Zb,0)
          return(power_result)
     }
     if(!missing(N1) & missing(N2)){
          Zb <- qnorm(power)
          N2_result <- (1-(R2*h2+c2)^2)^2 / (h2^2*(abs(R1-R2)^2)/((Za+Zb)^2) - ((1-(R1*h2+c2)^2)^2/N1))
          return(round(N2_result))
     }
     if(missing(N1) & !missing(N2)){
          Zb <- qnorm(power)
          N1_result <- (1-(R1*h2+c2)^2)^2 / (h2^2*(abs(R1-R2)^2)/((Za+Zb)^2) - ((1-(R2*h2+c2)^2)^2/N2))
          return(round(N1_result))
     }
     if(missing(N1) & missing(N2) & !is.null(p_N1)){
          Zb <- qnorm(power)
          N_total <- (1-(R1*h2+c2)^2)^2/(p_N1*h2^2*abs(R1-R2)^2/((Za+Zb)^2)) +  (1-(R2*h2+c2)^2)^2/((1-p_N1)*h2^2*(abs(R1-R2)^2)/((Za+Zb)^2))
          N1_result <- N_total * p_N1
          N2_result <- N_total * (1-p_N1)
          return(c(round(N1_result),round(N2_result)))

     }


}



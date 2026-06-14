#' kinsim_double2
#' @description The function to generate two groups of univariate kin pair(e.g., both MZ and DZ twins) data using a multivariate norm approach, given the ACE components.
#' \cr
#' \cr
#' Two approaches can be selected: a) simulate two groups of kin pairs using the genetic relatedness directly b) simulate two groups of kin pairs by combining MZ twins and DZ twins to achieve the required genetic relatedness (.5<R<1).
#' @param GroupNames A character vector specifying two names of the simulated kin pairs
#' @param GroupSizes A numeric vector specifying two group sizes indicating the amount of kin pairs in respective group.
#' @param GroupRel A numeric vector specifying two genetic relatedness values of the simulated kin pairs
#' @param GroupR_c A numeric vector specifying two common environment correlation coefficients of the simulated kin pairs
#' @param mu A numeric vector specifying two mean values for the generated variable of the kin pairs
#' @param ace1 A numeric vector specifying three variance components under an ACE (additive genetics, common environment, unique environment) structure for group1
#' @param ace2 A numeric vector specifying three variance components under an ACE (additive genetics, common environment, unique environment) structure for group2
#' @param prop_missing A numeric vector specifying the percentage random missing data for kin pairs
#' @param ifComb A logical value specifying the approach to achieve the required genetic relatedness value. \code{TRUE} = using combination approach. \code{FALSE} = using direct approach. (See function description for a detailed explanation of two approaches.)
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

kinsim_double2 <- function(GroupNames = c("KinPair1", "KinPair2"),
                          GroupSizes = c(100, 100),
                          GroupRel = c(1, .5),
                          GroupR_c = c(1, 1),
                          mu = c(0, 0),
                          ace1 = c(1, 1, 1),
                          ace2 = c(1, 1, 1),
                          prop_missing = c(.20,.20),
                          ifComb = FALSE) {
  if (ifComb==FALSE) {
    df_N1 <- kinsim_single(
      name = GroupNames[1],
      Rel = GroupRel[1],
      r_c = GroupR_c[1],
      n = GroupSizes[1],
      mu = mu[1],
      ace = ace1
    )
    df_N2 <- kinsim_single(
      name = GroupNames[2],
      Rel = GroupRel[2],
      r_c = GroupR_c[2],
      n = GroupSizes[2],
      mu = mu[2],
      ace = ace2
    )
    df_final <- rbind(df_N1, df_N2)
  } else {
    if ((GroupRel[1] == 1 | GroupRel[1] == .5) & GroupRel[2] != 1 & GroupRel[2] != .5) {

      print(paste("the if statement you think is running is running"))
      df_N1 <- kinsim_single(
        name = GroupNames[1],
        Rel = GroupRel[1],
        r_c = GroupR_c[1],
        n = GroupSizes[1],
        mu = mu[1],
        ace = ace1
      )
      df2MZ <- kinsim_single(
        name = GroupNames[2],
        Rel = 1,
        r_c = GroupR_c[2],
        n = round((GroupRel[2] - .5) * 2 * GroupSizes[2]),
        mu = mu[2],
        ace = ace2
      )
      df2DZ <- kinsim_single(
        name = GroupNames[2],
        Rel = .5,
        r_c = GroupR_c[2],
        n = GroupSizes[2] - round((GroupRel[2] - .5) * 2 * GroupSizes[2]),
        mu = mu[2],
        ace = ace2
      )
      df_N2 <- rbind(df2MZ, df2DZ)
      df_N2 <- df_N2[sample(1:nrow(df_N2)), ]
      df_N2$id <- 1:nrow(df_N2)
      df_N2$R <- GroupRel[2]
      df_final <- rbind(df_N1, df_N2)
    } else if (GroupRel[1] != 1 & GroupRel[1] != .5 & (GroupRel[2] == 1 | GroupRel[2] == .5)) {
      df1MZ <- kinsim_single(
        name = GroupNames[1],
        Rel = 1,
        r_c = GroupR_c[1],
        n = round((GroupRel[1] - .5) * 2 * GroupSizes[1]),
        mu = mu[1],
        ace = ace1
      )
      df1DZ <- kinsim_single(
        name = GroupNames[1],
        Rel = .5,
        r_c = GroupR_c[1],
        n = GroupSizes[1] - round((GroupRel[1] - .5) * 2 * GroupSizes[1]),
        mu = mu[1],
        ace = ace1
      )
      df_N1 <- rbind(df1MZ, df1DZ)
      df_N1 <- df_N1[sample(1:nrow(df_N1)), ]
      df_N1$id <- 1:nrow(df_N1)
      df_N1$R <- GroupRel[1]
      df_N2 <- kinsim_single(
        name = GroupNames[2],
        Rel = GroupRel[2],
        r_c = GroupR_c[2],
        n = GroupSizes[2],
        mu = mu[2],
        ace = ace2
      )
      df_final <- rbind(df_N1, df_N2)

    } else  if ((GroupRel[1] == 1 | GroupRel[1] == .5) & (GroupRel[2] == 1 | GroupRel[2] == .5)) {
      df_N1 <- kinsim_single(
        name = GroupNames[1],
        Rel = GroupRel[1],
        r_c = GroupR_c[1],
        n = GroupSizes[1],
        mu = mu[1],
        ace = ace1
      )
      df_N2 <- kinsim_single(
        name = GroupNames[2],
        Rel = GroupRel[2],
        r_c = GroupR_c[2],
        n = GroupSizes[2],
        mu = mu[2],
        ace = ace2
      )
      df_final <- rbind(df_N1, df_N2)

    } else if (GroupRel[1] != 1 & GroupRel[1] != .5 & GroupRel[2] != 1 & GroupRel[2] != .5) {

         print(paste("running the new condition"))

      df1MZ <- kinsim_single(
        name = GroupNames[1],
        Rel = 1,
        r_c = GroupR_c[1],
        n = round((GroupRel[1] - .5) * 2 * GroupSizes[1]),
        mu = mu[1],
        ace = ace1
      )
      df1DZ <- kinsim_single(
        name = GroupNames[1],
        Rel = .5,
        r_c = GroupR_c[1],
        n = GroupSizes[1] - round((GroupRel[1] - .5) * 2 * GroupSizes[1]),
        mu = mu[1],
        ace = ace1
      )
      df_N1 <- rbind(df1MZ, df1DZ)
      df_N1 <- df_N1[sample(1:nrow(df_N1)), ]
      df_N1$id <- 1:nrow(df_N1)
      df_N1$R <- GroupRel[1]

      df2MZ <- kinsim_single(
        name = GroupNames[2],
        Rel = 1,
        r_c = GroupR_c[2],
        n = round((GroupRel[2] - .5) * 2 * GroupSizes[2]),
        mu = mu[2],
        ace = ace2
      )
      df2DZ <- kinsim_single(
        name = GroupNames[2],
        Rel = .5,
        r_c = GroupR_c[2],
        n = GroupSizes[2] - round((GroupRel[2] - .5) * 2 * GroupSizes[2]),
        mu = mu[2],
        ace = ace2
      )
      df_N2 <- rbind(df2MZ, df2DZ)
      df_N2 <- df_N2[sample(1:nrow(df_N2)), ]
      df_N2$id <- 1:nrow(df_N2)
      df_N2$R <- GroupRel[2]

      df_final <- rbind(df_N1, df_N2)
    } else {
      stop("Invalid combination of GroupRel values. Please check the input.")
    }

  }
  # Add missingness and ordinal columns if prop_missing is greater than 0 for either group
  if (any(prop_missing > 0)) {
    df_final <- .add_missing_and_ordinal(df_final,
                                         GroupNames=GroupNames,
                                         prop_missing=prop_missing)

  }
  return(df_final) #I think the problem is here? my error message is about df_missing not being found - my only thought is that df_missing is the object returned by that function, but we only ask for df_final bc thats all we care about in the end
}

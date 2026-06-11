#' rmvn
#' @keywords internal
#' @param n Sample Size
#' @param sigma Covariance matrix
#' @return Generates multivariate normal data from a covariance matrix (\code{sigma}) of length \code{n}
#' @importFrom stats rnorm
rmvn <- function(n, sigma) {
  Sh <- with(
    svd(sigma),
    v %*% diag(sqrt(d)) %*% t(u)
  )
  matrix(stats::rnorm(ncol(sigma) * n),
    ncol = ncol(sigma)
  ) %*% Sh
}


# Apply random missingness by GroupName (not by R value, which breaks when groups
# share the same relatedness), then compute 4-category ordinal columns.
# Cutpoints: (-Inf, -2] = 1, (-2, -1] = 2, (-1, 1) = 3, [1, Inf) = 4.
# NA y values produce NA ordinal scores; missingness is cascaded from y1 to y2.
.add_missing_and_ordinal <- function(df, GroupNames, prop_missing) {

#' add missing
#' @description Internal function to add missingness to the simulated data frame.
#'
#' @details
#' Apply random missingness by GroupName (not by R value, which breaks when groups share the same relatedness),
#' @param df A data frame containing the simulated data with columns 'GroupName', 'y1', and 'y2'.
#' @param GroupNames A character vector of length 2 specifying the group names corresponding to the two groups in the data frame.
#' @param prop_missing A numeric vector of length 2 specifying the proportion of missing values
#' for each group. The first element corresponds to the first group in GroupNames, and the second element corresponds to the second group.
#' @return A modified data frame with missing values added to 'y1' and 'y2' according to the specified proportions, and new ordinal columns 'Ord_1' and 'Ord_2' added based on the cutpoints applied to 'y1' and 'y2 respectively.
#' @importFrom stats runif
.add_missing <- function(df, GroupNames, prop_missing) {

  in_g1 <- df$GroupName == GroupNames[1]
  in_g2 <- df$GroupName == GroupNames[2]

  miss_mask <- logical(nrow(df))
  miss_mask[in_g1] <- stats::runif(sum(in_g1)) < prop_missing[1]
  miss_mask[in_g2] <- stats::runif(sum(in_g2)) < prop_missing[2]

  df$y1[miss_mask] <- NA_real_
  df$y2[miss_mask] <- NA_real_

  df
}

#' add_missing_and_ordinal
#'
#' @details Applies random missingness to the 'y1' and 'y2' columns of the input data frame based on the specified proportions for each group, and
#' then compute 4-category ordinal columns.
#' Cutpoints: (-Inf, -2] = 1, (-2, -1] = 2, (-1, 1) = 3, [1, Inf) = 4.
#' NA y values produce NA ordinal scores; missingness is cascaded from y1 to y2.
#' @inheritParams .add_missing
#' @return A modified data frame with missing values added to 'y1' and 'y2' according to the specified proportions, and new ordinal columns 'Ord_1' and 'Ord_2' added based on the cutpoints applied to 'y1' and 'y2 respectively.
#'
.add_missing_and_ordinal <- function(df, GroupNames, prop_missing) {

   df_missing <-  .add_missing(df, GroupNames, prop_missing)


  .to_ord <- function(y) {
    out <- rep(NA_integer_, length(y))
    out[!is.na(y) & y <= -2]           <- 1L
    out[!is.na(y) & y > -2 & y <= -1]  <- 2L
    out[!is.na(y) & y > -1 & y <   1]  <- 3L
    out[!is.na(y) & y >=  1]           <- 4L
    out
  }


  df$Ord_1 <- .to_ord(df$y1)
  df$Ord_2 <- .to_ord(df$y2)
  df
}

  df_missing$Ord_1 <- .to_ord(df_missing$y1)
  df_missing$Ord_2 <- .to_ord(df_missing$y2)
  df_missing
}





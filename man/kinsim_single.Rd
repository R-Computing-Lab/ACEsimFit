% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinsim_single.R
\name{kinsim_single}
\alias{kinsim_single}
\title{kinsim_single}
\usage{
kinsim_single(
  name = "KinPair1",
  Rel = 1,
  r_c = 1,
  n = 100,
  mu = 0,
  ace = c(1, 1, 1)
)
}
\arguments{
\item{name}{Assigned name for the simulated group of kin pairs}

\item{Rel}{Genetic relatedness of the simulated kin pairs}

\item{r_c}{Assumed common environment correlation}

\item{n}{The number of generated kin pairs.(n PAIRS of data; The total number of participants is 2n)}

\item{mu}{The mean for generated variable}

\item{ace}{Vector of variance components under an ACE (additive genetics, common environment, unique environment) structure}
}
\value{
Returns \code{data.frame} with the following:
\item{GroupName}{group name of the kin pairs}
\item{R}{level of genetic relatedness for the kin pairs}
\item{r_c}{level of common envrionment correlation of the kin pairs}
\item{id}{id}
\item{A1}{Additive genetic component for kin1 of the kin pairs}
\item{A2}{Additive genetic component for kin2 of the kin pairs}
\item{C1}{shared-environmental component for kin1 of the kin pairs}
\item{C2}{shared-environmental component for kin2 of the kin pairs}
\item{E1}{non-shared-environmental component for kin1 of the kin pairs}
\item{E2}{non-shared-environmental component for kin2 of the kin pairs}
\item{y1}{generated variable i for kin1}
\item{y2}{generated variable i for kin2}
}
\description{
The function to generate one group of univariate kin pair (e.g., only DZ twins) data using a multivariate norm approach, given the ACE components.
}

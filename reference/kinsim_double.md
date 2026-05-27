# kinsim_double

The function to generate two groups of univariate kin pair(e.g., both MZ
and DZ twins) data using a multivariate norm approach, given the ACE
components.  
  
Two approaches can be selected: a) simulate two groups of kin pairs
using the genetic relatedness directly b) simulate two groups of kin
pairs by combining MZ twins and DZ twins to achieve the required genetic
relatedness (.5\<R\<1).

## Usage

``` r
kinsim_double(
  GroupNames = c("KinPair1", "KinPair2"),
  GroupSizes = c(100, 100),
  GroupRel = c(1, 0.5),
  GroupR_c = c(1, 1),
  mu = c(0, 0),
  ace1 = c(1, 1, 1),
  ace2 = c(1, 1, 1),
  ifComb = FALSE
)
```

## Arguments

- GroupNames:

  A character vector specifying two names of the simulated kin pairs

- GroupSizes:

  A numeric vector specifying two group sizes indicating the amount of
  kin pairs in respective group.

- GroupRel:

  A numeric vector specifying two genetic relatedness values of the
  simulated kin pairs

- GroupR_c:

  A numeric vector specifying two common environment correlation
  coefficients of the simulated kin pairs

- mu:

  A numeric vector specifying two mean values for the generated variable
  of the kin pairs

- ace1:

  A numeric vector specifying three variance components under an ACE
  (additive genetics, common environment, unique environment) structure
  for group1

- ace2:

  A numeric vector specifying three variance components under an ACE
  (additive genetics, common environment, unique environment) structure
  for group2

- ifComb:

  A logical value specifying the approach to achieve the required
  genetic relatedness value. `TRUE` = using combination approach.
  `FALSE` = using direct approach. (See function description for a
  detailed explanation of two approaches.)

## Value

Returns `data.frame` with the following:

- GroupName:

  group name of the kin pairs

- R:

  level of relatedness for the kin pair

- r_c:

  level of common envrionment correlation of the kin pairs

- id:

  id

- A1:

  Additive genetic component for kin1 of the kin pairs

- A2:

  Additive genetic component for kin2 of the kin pairs

- C1:

  shared-environmental component for kin1 of the kin pairs

- C2:

  shared-environmental component for kin2 of the kin pairs

- E1:

  non-shared-environmental component for kin1 of the kin pairs

- E2:

  non-shared-environmental component for kin2 of the kin pairs

- y1:

  generated variable i for kin1

- y2:

  generated variable i for kin2

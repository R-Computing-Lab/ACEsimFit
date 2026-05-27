# Sim_Fit

A function to simulate a set of kin pair data and fit them with ACE
models. Can be helpful with checking model performance for a given
parameter setting.

## Usage

``` r
Sim_Fit(
  GroupNames = c("KinPair1", "KinPair2"),
  GroupSizes = c(100, 100),
  nIter = 100,
  SSeed = 62,
  GroupRel = c(1, 0.5),
  GroupR_c = c(1, 1),
  mu = c(0, 0),
  ace1 = c(1, 1, 1),
  ace2 = c(1, 1, 1),
  ifComb = FALSE,
  lbound = FALSE,
  saveRaw = FALSE
)
```

## Arguments

- GroupNames:

  A character vector specifying two names of the simulated kin pairs

- GroupSizes:

  A numeric vector specifying two group sizes indicating the amount of
  kin pairs in respective group.

- nIter:

  A numeric value specifying the number of iteration you want to run
  given the parameters assigned (i.e. the number of model fitting
  results you want to get)

- SSeed:

  An integer specifying the starting seed of the random number. This
  parameter will make sure the simulated results are replicable across
  time

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

- lbound:

  A logical value indicating if a lower boundary of .0001 will be
  imposed to the estimated A, C and E components

- saveRaw:

  A logical value specifying if the raw simulated data should be saved
  in the output list

## Value

Returns a two-level `list`. Level-one is the number of iterations.
Level-two is the model fitting results and raw data (if
`saveRaw = TRUE`) of the simulated data from the respective iteration.
Level-two includes:

- Results:

  A `list` including 1) A `data.frame` displaying the nested comparison
  model between ACE, AE, CE, E models and 2) A `list` of all model fit
  information generated from OpenMx

- Data:

  A `data.frame` consists of the simulated raw data

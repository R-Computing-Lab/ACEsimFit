# kinsim_single

The function to generate one group of univariate kin pair (e.g., only DZ
twins) data using a multivariate norm approach, given the ACE
components.

## Usage

``` r
kinsim_single(
  name = "KinPair1",
  Rel = 1,
  r_c = 1,
  n = 100,
  mu = 0,
  ace = c(1, 1, 1)
)
```

## Arguments

- name:

  Assigned name for the simulated group of kin pairs

- Rel:

  Genetic relatedness of the simulated kin pairs

- r_c:

  Assumed common environment correlation

- n:

  The number of generated kin pairs.(n PAIRS of data; The total number
  of participants is 2n)

- mu:

  The mean for generated variable

- ace:

  Vector of variance components under an ACE (additive genetics, common
  environment, unique environment) structure

## Value

Returns `data.frame` with the following:

- GroupName:

  group name of the kin pairs

- R:

  level of genetic relatedness for the kin pairs

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

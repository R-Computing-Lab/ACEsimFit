# fit_uniACE

Use OpenMx to quickly fit a univariate ACE model

## Usage

``` r
fit_uniACE(
  data_1,
  data_2,
  GroupRel = c(1, 0.5),
  GroupR_c = c(1, 1),
  lbound = FALSE
)
```

## Arguments

- data_1:

  A n by 2 `data.frame` consisting of the group1 kin pairs

- data_2:

  A n by 2 `data.frame` consisting of the group2 kin pairs

- GroupRel:

  A numeric vector specifying two genetic relatedness values of two
  groups of kin pairs

- GroupR_c:

  A numeric vector specifying two common environment correlation
  coefficients of two groups of kin pairs

- lbound:

  A logical value indicating if a lower boundary of .0001 will be
  imposed to the estimated A, C and E components

## Value

Returns a `list` with the following:

- df_nested:

  A `data.frame` displaying the nested comparison model between ACE, AE,
  CE, E models

- fitACE:

  A `list` of all model fit information generated from OpenMx

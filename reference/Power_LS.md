# Power_LS

The function is designed for calculating the power of heritability
estimation from ACE models given the parameter settings. Or calculate
one of the parameter settings (N,R,h2,c2) given the rest of known
parameters. This power calculator is made based on the Least Squares
theory and follows the mathematical derivation proposed by Visscher
(2004).

## Usage

``` r
Power_LS(
  N1,
  N2,
  power,
  p_N1 = NULL,
  h2,
  c2,
  R1 = 1,
  R2 = 0.5,
  alpha = 0.05,
  digits = 0
)
```

## Arguments

- N1:

  The number of kin pairs for group1 (amount of PAIRS)

- N2:

  The number of kin pairs for group2

- power:

  The power of heritability estimation. Specified if you want to return
  the required sample sizes.

- p_N1:

  The proportion of kin group1 over the . Required to be specified if
  the user wants to calculate the N1 and N2 simultaneously.

- h2:

  The assumed standard heritability value of the target trait. 0 \< h2
  \< 1

- c2:

  The assumed standard common environmental effects on the target trait.
  0 \< c2 \< 2

- R1:

  The genetic relatedness of kin pair group1

- R2:

  The genetic relatedness of kin pair group2

- alpha:

  The type-one error rate for heritability estimation.

- digits:

  The number of digits to round the results. Default is 0.

## Value

A numeric `vector` of power when `N1` and `N2` are both specified.  
A numeric `vector` of `N1` (or `N2`) when `N2` (or `N1`) is specified. A
numeric `vector` of `N1` and `N2` when `RatioN` is specified.

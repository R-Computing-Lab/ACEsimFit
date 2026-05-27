# ACEsimFit

``` r

library(ACEsimFit)
```

## Overview

This package is designed for 1) Simulating kin pairs data based on the
assumption that every trait is affected by genetic effects (A), common
environmental effects (C) and unique environmental effects (E). 2) Using
kin pairs data to fit an ACE model and get model fit output. In our
vignette, we will use height as a trait to demonstrate how people can
simulate kin pairs of height data and how to use those simulated data to
fit a univariate ACE model. We will also discuss the possible ways to
analyse these model fitting results. 3) Calculate power of heritability
estimate given a specific condition. (Note: The ideas and data in this
vignette are only for demonstration purpose and can be inaccurate in
terms of the true facts in science.)

## Simulate kin pair data

Suppose we have a situation like this: I wants to investigate what
factors impact the height of the adults in a non-European country.
Researchers in developed countries like US and UK have collected a sea
of twin data to explore the question and reached a conclusion: 60% of
variance of height comes from genes, 20% of variance comes from family
environment and 20% comes from personal environment and error. I wanted
to see if this conclusion holds true in an Asian country. One approach
is to use a public dataset with height measure and family structure
records.

However, the dataset has two problems: 1) Only have 180 pairs of twins
2) can only distinguish MZ twins and same-sex DZ twins. Given the
problems this dataset had, it will be useful to do a priori power
analysis to check the statistical power of my study.

So first I need to simulate data which they have same variance and
covariance structure as the data in the public dataset. We can use the
`kinsim_double` function to do that. Here’s some conditions I want to
replicate:

- Sample sizes: 60 pairs of opposite-sex DZ twins (genetic relatedness
  is .5) and 120 pairs of same-sex MZ and DZ twins (genetic relatedness
  is .75)

- The relatedness of .75 is achieved by combining MZ and DZ twins
  instead of simulating kin pairs with relatedness of .75 directly.

- ACE variance structure: A = .6, C = .2, E = .2

- Group mean is both 0, because they should be Z scores.

Therefore, we set the parameters of the `kinsim_double` function like
below:

``` r

kindata <- kinsim_double(
  GroupNames = c("SStwins", "OStwins"),
  GroupSizes = c(120, 60),
  GroupRel = c(.75, 0.5),
  GroupR_c = c(1, 1),
  mu = c(0, 0),
  ace1 = c(.6, .2, .2),
  ace2 = c(.6, .2, .2),
  ifComb = TRUE
)
head(kindata)
#>    GroupName    R r_c id          A1          A2          C1          C2
#> 48   SStwins 0.75   1  1  0.02248168  0.02248168 -0.03448461 -0.03448461
#> 66   SStwins 0.75   1  2  0.87922301  0.73493224 -0.42473098 -0.42473098
#> 47   SStwins 0.75   1  3 -0.28043976 -0.28043976 -0.25064124 -0.25064124
#> 14   SStwins 0.75   1  4 -0.76368579 -0.76368579 -0.59577363 -0.59577363
#> 16   SStwins 0.75   1  5 -0.73345621 -0.73345621  0.18573124  0.18573124
#> 95   SStwins 0.75   1  6 -0.01048472 -0.06118696  0.19529037  0.19529037
#>             E1          E2          y1         y2
#> 48  0.09741996 -0.19555812  0.08541702 -0.2075610
#> 66 -0.49639316 -0.41603625 -0.04190114 -0.1058350
#> 47 -0.42105074 -0.03581002 -0.95213174 -0.5668910
#> 14 -0.33618089 -0.81363682 -1.69564031 -2.1730962
#> 16 -0.37056277  0.01884294 -0.91828774 -0.5288820
#> 95  0.39216748 -0.67556163  0.57697313 -0.5414582
```

Now you can see we have a data.frame with 180 pairs of simulated twins.
And we also have their variance components and other information.

## Calculate the power of the heritability estimate given the variance structure

Generally, calculating power for the A estimate had two approaches:
likelihood theory and the least squares theory.

### Calculate power based on likelihood theory

To use the likelihood theory, we need to simulate a number of datasets
with the same variance structure and average the suggested power from
each set of data. See a more detailed explanation at
[link](http://www.people.vcu.edu/~bverhulst/power/power.md). So we need
to have a set of model fitting results with the -2ll values for the ACE
model and the CE model.

Luckily, in our package they have a function `Sim_Fit` to simulate kin
pairs data, fit them automatically into a ACE model and return the model
summary results. The function fit the ACE model with the help of
`OpenMx` package.

Here, we again assigned the same parameters to the function. There are a
few new parameters here:

- We want 50 simulated datasets so my power calculation can be
  relatively accurate. So we set the `nIter = 50`. Set more iterations
  when you want to get a even more accurate power.
- We hope my simulated results can be replicated so we set `SSeed = 62`.
  Pick a lucky number here!
- We don’t want to contrain the estimation so we set `lbound = FALSE`.
- We don’t need the raw data since they will eat up my poor-RAM laptop.
  So we set `saveRaw = FALSE`.

``` r

time1 <- Sys.time()
results_fit <- Sim_Fit(
  GroupNames = c("SStwins", "OStwins"),
  GroupSizes = c(120, 60),
  nIter = 50,
  SSeed = 62,
  GroupRel = c(.75, 0.5),
  GroupR_c = c(1, 1),
  mu = c(0, 0),
  ace1 = c(.6, .2, .2),
  ace2 = c(.6, .2, .2),
  ifComb = TRUE,
  lbound = FALSE,
  saveRaw = FALSE
)
time2 <- Sys.time()
## FYI, the time used for the results above is here. So design your simulation wisely!!!
time2 - time1
#> Time difference of 43.02106 secs
```

Here’s one example of the nested comparison table from the results

``` r

results_fit[["Iteration1"]][["Results"]][["nest"]]
#>            base comparison ep minus2LL  df      AIC     diffLL diffdf
#> 1 oneACEvc_1cov       <NA>  4 896.1042 356 904.1042         NA     NA
#> 2 oneACEvc_1cov    oneAEvc  3 897.0780 357 903.0780  0.9737667      1
#> 3 oneACEvc_1cov    oneCEvc  3 898.0232 357 904.0232  1.9190034      1
#> 4 oneACEvc_1cov     oneEvc  2 980.2935 358 984.2935 84.1892959      2
#>              p
#> 1           NA
#> 2 3.237426e-01
#> 3 1.659666e-01
#> 4 5.230301e-19
```

We can then calculated the weighted ncp with the average difference of
log-likelihood for ACE and CE models.In turn, I calculated the power for
the given variance structure and relatedness:

``` r

N <- 180 ## the total number of kin pairs you used in your previous simulation
## Calculate the average diffLL between ACE and CE model.
DiffLL <- numeric()
for (i in 1:50) {
  DiffLL[i] <- results_fit[[1]][["Results"]][["nest"]]$diffLL[3]
}
meanDiffLL <- mean(DiffLL)
## Calculate the power based on an alpha level of .05
Power <- 1 - pchisq(qchisq(1 - .05, 1), 1, meanDiffLL)
Power
#> [1] 0.2831639
```

So you can see the power for my model would be insufficient to have a
confident estimate of height heritability. :(

### Calculate power based on least squares theory

Unlike the likelihood theory, we don’t need to run simulations to
calculate power because we can have a formula in Xuanyu and Mason’s
paper (submitted; you will see a link here soon). The function
`Power_LS` is designed for power calculation based on LS theory.

Here I typed the parameters my research question involved into the
functions.

``` r

Power_LS(N1 = 120, N2 = 60, h2 = .6, c2 = .2, R1 = .75, R2 = 0.5, alpha = 0.05)
#> [1] 0.3881047
```

Your can see the power calculated with LS theory is different from the
one with likelihood theory. This is due to 1) with the LS formula, we
did’t consider the combination of MZ and DZ twins to reach a relatedness
of .75 2) they will be different especially when the sample size is
small. In our case, if we have 1800 pairs in stead of 180, the power
from likelihood theory will be 0.9890003 and the power from least
squares theory will be 0.9960664. (Much more similar at this time!)

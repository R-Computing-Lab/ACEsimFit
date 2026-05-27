# Changelog

## ACEsimFit (development version)

- Added a `NEWS.md` file to track changes to the package.
- Exposed rounding parameter in
  [`Power_LS()`](https://r-computing-lab.github.io/ACEsimFit/reference/Power_LS.md)
  via the `digits` argument.
- Fixed silent `NULL` return in
  [`Power_LS()`](https://r-computing-lab.github.io/ACEsimFit/reference/Power_LS.md)
  when both `N1` and `N2` are missing and `p_N1` is unspecified;
  function now stops with an informative error message.
- Fixed parameter shadowing in
  [`kinsim_single()`](https://r-computing-lab.github.io/ACEsimFit/reference/kinsim_single.md)
  where the `r_c` input was overwritten by an internal vector of the
  same name.
- Fixed `id` column generation in
  [`kinsim_single()`](https://r-computing-lab.github.io/ACEsimFit/reference/kinsim_single.md)
  to use [`seq_len()`](https://rdrr.io/r/base/seq.html), avoiding the
  `1:0` pitfall with edge-case inputs.
- Replaced `T` with `TRUE` in
  [`fit_uniACE()`](https://r-computing-lab.github.io/ACEsimFit/reference/fit_uniACE.md)
  sub-model runs to guard against accidental overwriting of `T`.
- Expanded test suite: added structural, correctness, and statistical
  tests for
  [`kinsim_single()`](https://r-computing-lab.github.io/ACEsimFit/reference/kinsim_single.md),
  [`kinsim_double()`](https://r-computing-lab.github.io/ACEsimFit/reference/kinsim_double.md),
  [`Power_LS()`](https://r-computing-lab.github.io/ACEsimFit/reference/Power_LS.md),
  and
  [`Sim_Fit()`](https://r-computing-lab.github.io/ACEsimFit/reference/Sim_Fit.md).

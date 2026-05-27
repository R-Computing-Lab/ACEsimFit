# ACEsimFit (development version)
* Added a `NEWS.md` file to track changes to the package.
* Exposed rounding parameter in `Power_LS()` via the `digits` argument.
* Fixed silent `NULL` return in `Power_LS()` when both `N1` and `N2` are missing and `p_N1` is unspecified; function now stops with an informative error message.
* Fixed parameter shadowing in `kinsim_single()` where the `r_c` input was overwritten by an internal vector of the same name.
* Fixed `id` column generation in `kinsim_single()` to use `seq_len()`, avoiding the `1:0` pitfall with edge-case inputs.
* Replaced `T` with `TRUE` in `fit_uniACE()` sub-model runs to guard against accidental overwriting of `T`.
* Expanded test suite: added structural, correctness, and statistical tests for `kinsim_single()`, `kinsim_double()`, `Power_LS()`, and `Sim_Fit()`.

# ACEsimFit (0.0.0.9)
* Uploaded to CRAN on 2022-10-10 22:03:05 UTC

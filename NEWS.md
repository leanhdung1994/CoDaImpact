# CoDaImpact (development version)

## Bug fixes

* `VariationTable()` computed the new value of a scalar response as
  `y0 * (1 + h * eta)` instead of `y0 + h * eta`, so its rows "New value",
  "Variation in %" and "Variation in units" were wrong for a scalar response.
* `VariationTable()` now stops with its own message when neither `Xvar` nor
  the response is compositional, and rejects invalid `obs`, `inc_size`,
  `inc_rate`, `Ytotal` and `normalize_Xdir`. It warns when `inc_rate` is
  given with a numeric `Xdir`, where it is ignored.

## Breaking changes

* `VariationScenario()` always names the columns of its result `Y` and `X`;
  a scalar variable kept its own name before (e.g. `VS$YIELD`). The names of
  the variables are in the attributes `Yvar` and `Xvar`, and a print method
  shows them. A scalar response is now a vector instead of a one-column
  matrix.

## Bug fixes

* `VariationScenario()`: `normalize_Xdir = FALSE` had no effect for a
  compositional `Xvar`; the direction was always scaled to norm 1.
* `VariationScenario()` rejects invalid `obs`, `inc_size`, `n_steps` and
  `normalize_Xdir`, and stops with a clear message when `Xdir` is missing for
  a compositional `Xvar`.

## Minor improvements and fixes

* A numeric direction `Xdir` (in `ShareRatioElasticities()`,
  `VariationScenario()` and `VariationTable()`) must have positive parts that
  are not all equal; zeros, `NA` and the neutral element are rejected with a
  clear error.

## Bug fixes

* The clr covariances `VARCOV_CLR` computed by `transformationSummary()`, and
  therefore the standard errors and intervals of `confint.lmCoDa()`, were
  wrong for variables in alr coordinates. Results for ilr coordinates do not
  change.

## Minor improvements and fixes

* `ShareRatioElasticities()` matches a named numeric `Xdir` to the parts of
  `Xvar` by name, and stops when the names are not the parts of `Xvar`
  (they were matched by position).

## Breaking changes

* `Impacts()` is renamed `Elasticities()`. The new `Impacts()` computes
  pairwise share ratio impacts (equations (29), (32) and (35) of Dargel and
  Thomas-Agnan, 2024) and has the arguments `object`, `Xvar` and `Xdir`.
  Code that calls `Impacts(object, Xvar, obs)` must use `Elasticities()`.

## Minor improvements and fixes

* `confint.lmCoDa()` rejects a non-numeric `level` and a `y_ref` that is not a
  single number or name, and warns when `y_ref` or `obs` is ignored because
  the response is scalar.
* Typos fixed in messages of `confint.lmCoDa()`; the attribute `type` of its
  result now reads "CI for clr coefficients".

## Documentation

* The help pages of `rice_yields`, `election`, `car_market` and
  `toulouse_retail` say what one row is and describe all columns; `election`
  now documents `unemp_rate`.

## Breaking changes

* `CoDa_seq()` and `CoDa_path()` return compositions closed to 1; rows are no
  longer rescaled to the total of `comp_from`.
* `CoDa_seq()` and `CoDa_path()` replace zeros by `1e-10` after closing, so
  results for inputs with zeros change.
* `CoDa_seq()` stops when `comp_from` and `comp_to` are the same composition.

## Bug fixes

* `CoDa_seq()` and `CoDa_path()` reject invalid inputs (`NA`, negative parts,
  non-logical `add_opposite`, non-positive `step_size`, ...) with clear errors.
* The default end point of `CoDa_seq()` no longer depends on `comp_from`.

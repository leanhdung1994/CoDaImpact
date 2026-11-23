# CoDaImpact (development version)

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

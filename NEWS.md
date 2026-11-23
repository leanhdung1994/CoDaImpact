# CoDaImpact (development version)

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

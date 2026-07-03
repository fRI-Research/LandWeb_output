Known issues: <https://github.com/fRI-Research/LandWeb_output/issues>

# LandWeb_output 2.0.0 (2023-09-21)

Note: this module is being retired in LandWeb v3, superseded by the
`NRV_summary` module together with the `nrvtools` package. This is the final
release of `LandWeb_output`.

* migrated from `raster` to `terra`: all `expectsInput`/`createsOutput` layers
  changed from `Raster`/`RasterLayer`/`RasterStack` to `SpatRaster`, added
  `terra` to `reqdPkgs`, and switched the prepInputs loader from
  `raster::raster` to `terra::rast`.
* new `standAgeMap` output: stand ages are now derived from `cohortData` via
  `standAgeMapGenerator()` (previously only `vegTypeMap` was produced).
* stand ages are computed biomass-weighted (`weight = "biomass"`), and this is
  now made explicit.
* added output descriptions (`standAgeMap`, `vegTypeMap`).
* bumped dependency requirements: `LandR` (>= 1.1.0.9063) and
  `SpaDES.tools` (>= 2.0.0).
* added `loadOrder` metadata (run after `Biomass_regeneration` and
  `Biomass_regenerationPM`).
* metadata fixes (package branch typo), added a `render-module-rmd` GitHub
  Action, and misc cleanup.

# LandWeb_output 1.3.3 (2022-09-19)

* plotting fix.
* trimmed the version metadata list down to the module itself (dropped pinned
  `LandR`/`SpaDES.core` entries).

# LandWeb_output 1.3.2 (2018-12-24)

* set a minimum `SpaDES.core` version.

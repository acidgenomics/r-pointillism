## pointillism 0.4.7 (2020-01-30)

### Minor changes

- Now using cli package for messages.
- Updated basejump dependencies.

## pointillism 0.4.6 (2020-01-20)

### Minor changes

- Reworked NAMESPACE to support migration of bioverbs package to acidgenerics.

## pointillism 0.4.5 (2020-01-03)

### Major changes

- `normalize`: reworked internal code to reflect changes in scater package.
  Instead of calling `normalizeSCE` or `normalize` internally, `normalizeCounts`
  is the recommended method for returning `normcounts` and `logcounts` assays
  that can be slotted in `assays` of the `SingleCellExperiment` object.
- Deprecated `plotDot` in favor of `plotDots` (note plural spelling change).
  This reads better and improves consistency with conventions used in scater.

### Minor changes

- Removed user-defined option to reduce verbosity, where applicable. If these
  functions are too noisy, we'll tone it down in the future.

## pointillism 0.4.4 (2019-10-30)

Bioconductor 3.10 compatibility fixes.

## pointillism 0.4.3 (2019-10-26)

### Minor changes

- Bug fix for legacy Seurat objects created via pointillism with named rownames
  containing gene identifiers mapped to the gene symbols. Now the coercion
  method for `Seurat` to `SingleCellExperiment` checks for this and sanitizes
  the names automatically so the resulting `SingleCellExperiment` is valid.

## pointillism 0.4.2 (2019-09-04)

### Major changes

- Reworked internal code to remove dplyr, magrittr, rlang, tibble, and tidyr
  dependencies. Consistently using base R and Bioconductor S4 methods internally
  instead, working primarily with `DataFrame` instead of `tbl_df`.

## pointillism 0.4.1 (2019-08-12)

### Minor changes

- Bug fix for gene symbols being sanitized by `as.data.frame` or `as_tibble`
  coercion inconsitently across R versions. Switched to using `as` coercion
  from `DataFrame` to `data.frame` instead.
- Improved message consistency, removing use of backticks in favor of single
  quotes, which is now the convention across basejump.
- Lintr check fixes.

## pointillism 0.4.0 (2019-08-07)

Initial support for monocle3 `cell_data_set` class.

### Major changes

- Reworked internal S4 function names.
- Updated internal example data.
- Moved cell-cycle and cell-type markers from Google Sheets here into the
  package internally at `inst/extdata/markers/`.
- Working on preliminary support for monocle3 `cell_data_set`.

## pointillism 0.3.2 (2019-04-25)

### Minor changes

- S4 generic reexport documentation fixes.

## pointillism 0.3.1 (2019-04-23)

### Minor changes

- Switch to importing ggplot2 code from new [acidplots][] package. Recommended
  default ggplot2 themes are now named `acid_theme_light` and `acid_theme_dark`.
- Consolidated S4 class and generator function documentation into single Rd
  files, where applicable.

## pointillism 0.3.0 (2019-04-16)

### Major changes

- Updated S4 methods to support [Seurat][] v3 release. Also updated imports
  from Seurat to reflect v3 change: `GetAssayData`, `Stdev`, `VariableFeatures`.
- Now importing S3 coercion methods from Seurat: `as.SingleCellExperiment`.
  `as.Seurat`, which coerces `SingleCellExperiment`, isn't particularly helpful
  because it requires `logcounts` slot to be defined in SCE. We're calling
  `CreateSeuratObject` internally in our S4 coercion method instead, so that
  any SCE object can be easily coerced to a `Seurat` object.
- Resaved and renamed example data: `seurat`, `seurat_all_markers`,
  `seurat_known_markers`.
- Overhauled R Markdown templates.

### Minor changes

- Added preliminary `plotCellTypesPerCluster` method support. Need to test this
  out on larger datasets with many samples per cluster.
- Switched Travis CI configuration to use Docker image.
- Added some additional comments regarding reticulate umap-learn configuration
  required to run `Seurat::RunUMAP`. Consider putting some information in the
  README, since there is a lot of misinformation on GitHub regarding proper
  UMAP configuration.
- Removed `setOldClass` definitions for `session_info` and `package_version`.
- Renamed `as-methods.R` to `coerce-methods.R`, matching bsaejump conventions.
- `cellCountsPerCluster`: Improved internal `interestingGroups` handling.
- Reworked reexport method for S4 generics imported from [bioverbs][] package.

## pointillism 0.2.5 (2019-04-01)

### Major changes

- Migrated code to [Acid Genomics][].
- Made `runZinbwave` defunct. No longer recommending zinbwave calculations for
  droplet single-cell RNA-seq.
- ggplot2 themes `theme_paperwhite` and `theme_midnight` are now imported from
  minimalism package instead of basejump.

## pointillism 0.2.4 (2019-01-24)

### Minor changes

- `SeuratMarkersPerCluster` to `tbl_df` coerce needed a `decode` call to
  properly handling gene-to-symbol mappings with run-length encoding (`Rle`).

## pointillism 0.2.3 (2019-01-08)

### Minor changes

- Comment out lines to pass lintr checks.
- Split out imports into a separate `imports.R` file.
- Added `nullOK` flag in assert checks, where applicable. See `isGGScale`, for
  example.

## pointillism 0.2.2 (2018-12-22)

### Major changes

- Migrated S4 generics to [bioverbs][] package.
- Finish migrating assert checks to [goalie][] package.

### Minor changes

- Moved S4 validity checks into `setClass` call, rather than using a separate
  `setValidity` call. Refer to `AllClasses.R` for details.
- Removed `makeCellTypeMarkersFromGoogle` utility function.

## pointillism 0.2.1 (2018-12-01)

### Major changes

- Reworked internal code to use [goalie][] package for assert checks.
- `SeuratMarkersPerCluster`: Simplified S4 validity checks.
- Improved `show` methods for S4 classes: `CellCycleMarkers`, `CellTypeMarkers`.
- Added `summary` method support for `SeuratMarkers`, `SeuratMarkersPerCluster`.

### Minor changes

- Miscellaneous documentation improvements and fixes.
- Reworked internal code shared between `CellCycleMarkers` and
  `CellTypeMarkers` S4 generator functions.
- `SeuratMarkers`: Overhauled internal code, moving from `.seuratMarkers`
  internal function. Refer to `SeuratMarkers.R` for details.
- Added AppVeyor CI support.

## pointillism 0.2.0 (2018-11-15)

### New functions

- `CellCycleMarkers`, `CellTypeMarkers`, and `KnownMarkers` generator functions.
- `cellTypesPerCluster` S4 generic.
- `topMarkers` S4 generic.
- New functions to import sheets from Google: `importCellCycleMarkersFromGoogle`
  and `importCellTypeMarkersFromGoogle`. These may be removed in a future
  update, in favor of simply using CSV files managed in the package. We've also
  included a helpful `makeCellTypeMarkersFromGoogle` utility function.
- `SeuratMarkers`: New function for returning sanitized Seurat marker output.
- `SeuratMarkersPerCluster`: Returns an S4 class split by cluster.

### Major changes

- Improved S4 coercion methods. Added support for `CellCycleMarkers` to
  `tbl_df`, `CellTypeMarkers` to `tbl_df`, and `SeuratMarkersPerCluster` to
  `tbl_df`.
- Reworked internal code for `plotFeature`.
- `plotKnownMarkers` is now exported instead of `plotKnownMarkersDetected`.
- Improved formals (and consistency) for dimensional reduction plotting
  functions. This applies to `plotReducedDim`, `plotTSNE`, `plotPCA`, and
  `plotUMAP`.
- `runZinbwave`: Improved internal code and weight handling.

### Minor changes

- Reworked validity checks for S4 classes.
- Reworked internal assert checks. Refer to `.assertIsBPPARAM`, for example.
- Broke out S4 methods into internal functions, where applicable. Refer to
  `cellCountsPerCluster.SingleCellExperiment`, for example.
- Cleaned up example datasets. Removed `known_markers_detected_small` and
  `sce_small`.
- `diffExp`: Reorganized and reworked internal code that uses zinbwave.
- Reworked and improved Seurat-to-SingleCell interconversion support. Refer
  to `seurat-SingleCellExperiment-methods` for details.
  
### Deprecations

These functions have been deprecated: `knownMarkers`, `knownMarkersDetected`,
`plotKnownMarkersDetected`, `readCellTypeMarkers`, `sanitizeMarkers`.

## pointilism 0.1.3 (2018-09-19)

### Major changes

- `as(seurat, "SingleCellExperiment")` now returns `scaleData` in assays, along
  with `varGenes` in `metadata` and ensures that `colData` and `rowRanges`
  are sanitized into camel case using the `camel` function internally.

### Minor changes

- Removed `convertGenesToSymbols` support for `seurat` class. Coerce genes
  to symbols manually prior to `seurat` coercion from `SingleCellExperiment`.
- `clusterID`: improved consistency between `SingleCellExperiment` and
  `seurat` methods, which now both work on `colData`.
- Improved internal `fetchData` calls to use camel case for consistency.
- `plotFeature` now sanitizes the `features` argument into camel case
  internally, for more consistent matching.
- `plotReducedDim` now uses camel case formatting for axis labels.
- Removed unnecessary legacy `.minimalAxis` constructor.
- `metadata<-` assignment method now stashes directly into `object@misc` for
  `seurat` class.
- Updated example data script and resaved.

## pointillism 0.1.2 (2018-08-22)

### New functions

- `clusterID`: S4 generic that returns the cell cluster identity mappings as a
  `factor` for all cells.
- `findMarkers`: Utility function that wraps `diffExp` to identify cluster-
  specific genes across all clusters.

### Minor changes

- Converted `runZinbwave` to an S4 method that works on
  `SingleCellExperiment`. For `seurat` objects, coerce to `SingleCellExperiment`
  before calculating weights.
- `plotDot` now uses default [ggplot2][] color scale, unless specified.
- `plotReducedDim`, `plotMarker`: `theme_paperwhite` is no longer used by
  default. Removed `grid` and `aspectRatio` arguments to simplify [ggplot2][]
  code handling.

## pointillism 0.1.1 (2018-08-20)

### New functions

- Now exporting `runZinbwave` and `hasZinbwave`, which were previously internal
  functions called by `diffExp`.

### Minor changes

- Improved internal handling of `assays` for objects that have zinbwave
  calculations applied. Now it doesn't attempt to keep only the `counts` assay.

## pointillism 0.1.0 (2018-08-10)

- Initial release.

[acidplots]: https://acidplots.acidgenomics.com/
[Acid Genomics]: https://acidgenomics.com/
[basejump]: https://basejump.acidgenomics.com/
[bioverbs]: https://bioverbs.acidgenomics.com/
[ggplot2]: https://ggplot2.tidyverse.org/
[goalie]: https://goalie.acidgenomics.com/
[minimalism]: https://minimalism.acidgenomics.com/
[Seurat]: https://satijalab.org/seurat/

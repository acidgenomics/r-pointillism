## """
## Cell-cycle and cell-type markers.
## Updated 2020-10-12.
##
## Here we're matching the stored Ensembl identifiers (`geneID`) using
## ensembldb to obtain the latest symbol names from Ensembl.
##
## Cell-cycle markers devied from:
## - Tirosh et al, 2015
## - http://satijalab.org/seurat/cell_cycle_vignette.html
## """

library(usethis)
library(stringr)

markers_dir <- system.file(
    "extdata", "markers",
    package = "pointillism",
    mustWork = TRUE
)
stopifnot(dir.exists(markers_dir))

## Ensembl release version.
release <- as.integer(readLines(file.path(markers_dir, "ensembl-release.txt")))

## Cell-cycle markers ==========================================================
cell_cycle_dir <- file.path(markers_dir, "cell-cycle")
stopifnot(dir.exists(cell_cycle_dir))
files <- list.files(path = cell_cycle_dir, pattern = "*.csv", full.names = TRUE)
cell_cycle_markers_list <- lapply(
    X = files,
    FUN = function(file) {
        importCellCycleMarkers(
            file = file,
            organism = basenameSansExt(file),
            release = release
        )
    }
)
names <- camelCase(basenameSansExt(files))
names(cell_cycle_markers_list) <- names

## Cell-type markers ===========================================================
cell_type_dir <- file.path(markers_dir, "cell-type")
stopifnot(dir.exists(cell_type_dir))
files <- list.files(path = cell_type_dir, pattern = "*.csv", full.names = TRUE)
cell_type_markers_list <- lapply(
    X = files,
    FUN = function(file) {
        importCellTypeMarkers(
            file = file,
            organism = basenameSansExt(file),
            release = release
        )
    }
)
names <- camelCase(basenameSansExt(files))
names(cell_type_markers_list) <- names

## Save R data ==================================================================
use_data(
    cell_cycle_markers_list,
    cell_type_markers_list,
    compress = "xz",
    overwrite = TRUE
)

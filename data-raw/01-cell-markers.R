## """
## Cell-cycle and cell-type markers.
## Updated 2021-03-03.
##
## Here we're matching the stored Ensembl identifiers (`geneId`) using
## ensembldb to obtain the latest symbol names from Ensembl.
##
## Cell-cycle markers devied from:
## - Tirosh et al, 2015
## - http://satijalab.org/seurat/cell_cycle_vignette.html
## """

library(usethis)
library(stringr)

markersDir <- system.file(
    "extdata", "markers",
    package = "pointillism",
    mustWork = TRUE
)
stopifnot(dir.exists(markersDir))

## Ensembl release version.
release <- as.integer(readLines(file.path(markersDir, "ensembl-release.txt")))

## Cell-cycle markers ==========================================================
cellCycleDir <- file.path(markersDir, "cell-cycle")
stopifnot(dir.exists(cellCycleDir))
files <- sort(list.files(
    path = cellCycleDir,
    pattern = "*.csv",
    full.names = TRUE
))
cellCycleMarkersList <- lapply(
    X = files,
    FUN = function(file) {
        organism <- sentenceCase(gsub("-", " ", basenameSansExt(file)))
        importCellCycleMarkers(
            file = file,
            organism = organism,
            release = release,
            ignoreVersion = TRUE
        )
    }
)
names <- camelCase(basenameSansExt(files))
names(cellCycleMarkersList) <- names

## Cell-type markers ===========================================================
cellTypeDir <- file.path(markersDir, "cell-type")
stopifnot(dir.exists(cellTypeDir))
files <- sort(list.files(
    path = cellTypeDir,
    pattern = "*.csv",
    full.names = TRUE
))
cellTypeMarkersList <- lapply(
    X = files,
    FUN = function(file) {
        organism <- sentenceCase(gsub("-", " ", basenameSansExt(file)))
        importCellTypeMarkers(
            file = file,
            organism = organism,
            release = release,
            ignoreVersion = TRUE
        )
    }
)
names <- camelCase(basenameSansExt(files))
names(cellTypeMarkersList) <- names

## Save R data ==================================================================
use_data(
    cellCycleMarkersList,
    cellTypeMarkersList,
    compress = "xz",
    overwrite = TRUE
)

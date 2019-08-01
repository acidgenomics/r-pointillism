## Cell-cycle and cell-type markers.
## Updated 2019-07-31.

## Here we're matching the stored Ensembl identifiers (`geneID`) using
## ensembldb to obtain the latest symbol names from Ensembl.

## Cell-cycle markers devied from:
## - Tirosh et al, 2015
## - http://satijalab.org/seurat/cell_cycle_vignette.html

library(usethis)
library(stringr)

# Import as `DataFrame` instead of `tbl_df`.
options("acid.data.frame") <- "DataFrame"

markersDir <- system.file(
    file.path("inst", "extdata", "markers"),
    package = "pointillism"
)
stopifnot(dir.exists(markersDir))

## Ensembl release version.
release <- as.integer(readLines(file.path(markersDir, "ensembl-release.txt")))

## Cell-cycle markers ==========================================================
cellCycleDir <- file.path(markersDir, "cell-cycle")
stopifnot(dir.exists(cellCycleDir))
files <- list.files(path = cellCycleDir, pattern = "*.csv", full.names = TRUE)
cellCycleMarkersList <- lapply(
    X = files,
    FUN = function(file) {
        data <- import(file = file)
        data <- as(data, "DataFrame")
        organism <- basenameSansExt(file)
        gene2symbol <- makeGene2SymbolFromEnsembl(
            organism = organism,
            release = release
        )
        CellCycleMarkers(
            object = data,
            gene2symbol = gene2symbol
        )
    }
)
names <- camel(basenameSansExt(files))
names(cellCycleMarkersList) <- names

## Cell-type markers ===========================================================
cellTypeDir <- file.path(markersDir, "cell-type")
stopifnot(dir.exists(cellTypeDir))
files <- list.files(path = cellTypeDir, pattern = "*.csv", full.names = TRUE)
cellTypeMarkersList <- lapply(
    X = files,
    FUN = function(file) {
        data <- import(file = file)
        data <- as(data, "DataFrame")
        organism <- basenameSansExt(file)
        gene2symbol <- makeGene2SymbolFromEnsembl(
            organism = organism,
            release = release
        )
        CellTypeMarkers(
            object = data,
            gene2symbol = gene2symbol
        )
    }
)
names <- camel(basenameSansExt(files))
names(cellTypeMarkersList) <- names

## Save R data ==================================================================
use_data(
    cellCycleMarkersList,
    cellTypeMarkersList,
    compress = "xz",
    overwrite = TRUE
)

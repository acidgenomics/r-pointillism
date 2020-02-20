#' Import cell-cycle markers
#'
#' @export
#' @note Updated 2020-02-20.
#'
#' @inheritParams acidroxygen::params
#' @inheritParams basejump::makeGene2SymbolFromEnsembl
#'
#' @examples
#' markers_dir <- system.file(
#'     file.path("inst", "extdata", "markers"),
#'     package = "pointillism"
#' )
#'
#' cell_cycle_dir <- file.path(markers_dir, "cell-cycle")
#' stopifnot(dir.exists(cell_cycle_dir))
#' files <- list.files(path = cell_cycle_dir, pattern = "*.csv", full.names = TRUE)
#' file <- files[[1L]]
#'
#' ## Ensembl release version.
#' release_file <- file.path(markers_dir, "ensembl-release.txt")
#' release <- as.integer(readLines(release_file))
#'
#' importCellCycleMarkers(
#'     file = file,
#'     organism = basenameSansExt(file),
#'     release = release
#' )
importCellCycleMarkers <- function(
    file,
    organism,
    release
) {
    object <- import(file)
    assert(isSubset(c("phase", "geneID", "modified"), colnames(object)))
    object <- as(object, "DataFrame")
    gene2symbol <- makeGene2SymbolFromEnsembl(
        organism = organism,
        release = release
    )
    CellCycleMarkers(
        object = object,
        gene2symbol = gene2symbol
    )
}

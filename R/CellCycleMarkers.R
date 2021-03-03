#' Cell-cycle markers
#'
#' @name CellCycleMarkers
#' @note Updated 2021-03-02.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams AcidGenomes::makeGene2SymbolFromEnsembl
#'
#' @return `CellCycleMarkers`.
#'
#' @examples
#' markersDir <- system.file(
#'     file.path("extdata", "markers"),
#'     package = "pointillism"
#' )
#'
#' cellCycleDir <- file.path(markersDir, "cell-cycle")
#' files <- list.files(cellCycleDir, pattern = "*.csv", full.names = TRUE)
#' file <- files[[1L]]
#'
#' organism <- sentenceCase(gsub("-", " ", basenameSansExt(file)))
#'
#' ## Ensembl release version.
#' releaseFile <- file.path(markersDir, "ensembl-release.txt")
#' release <- as.integer(readLines(releaseFile))
#'
#' importCellCycleMarkers(
#'     file = file,
#'     organism = organism,
#'     release = release
#' )
NULL



#' @rdname CellCycleMarkers
#' @export
CellCycleMarkers <-  # nolint
    function(object, gene2symbol) {
        assert(is(object, "DataFrame"))
        class <- "CellCycleMarkers"
        data <- .CellMarkers(
            object = object,
            gene2symbol = gene2symbol,
            class = class
        )
        new(Class = class, data)
    }



#' @rdname CellCycleMarkers
#' @export
importCellCycleMarkers <- function(
    file,
    organism,
    release
) {
    object <- import(file)
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

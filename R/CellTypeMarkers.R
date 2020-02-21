#' Cell-type markers
#'
#' @name CellTypeMarkers
#' @note Updated 2020-02-21.
#'
#' @inheritParams acidroxygen::params
#' @inheritParams basejump::makeGene2SymbolFromEnsembl
#'
#' @return `CellTypeMarkers`.
#'
#' @examples
#' markers_dir <- system.file(
#'     file.path("extdata", "markers"),
#'     package = "pointillism"
#' )
#'
#' cell_type_dir <- file.path(markers_dir, "cell-type")
#' files <- list.files(cell_type_dir, pattern = "*.csv", full.names = TRUE)
#' file <- files[[1L]]
#'
#' ## Ensembl release version.
#' release_file <- file.path(markers_dir, "ensembl-release.txt")
#' release <- as.integer(readLines(release_file))
#'
#' importCellTypeMarkers(
#'     file = file,
#'    organism = basenameSansExt(file),
#'     release = release
#' )
NULL



#' @rdname CellTypeMarkers
#' @export
CellTypeMarkers <-  # nolint
    function(object, gene2symbol) {
        assert(is(object, "DataFrame"))
        class <- "CellTypeMarkers"
        data <- .CellMarkers(
            object = object,
            gene2symbol = gene2symbol,
            class = class
        )
        new(Class = class, data)
    }



#' @rdname CellTypeMarkers
#' @export
importCellTypeMarkers <- function(
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
    CellTypeMarkers(
        object = object,
        gene2symbol = gene2symbol
    )
}

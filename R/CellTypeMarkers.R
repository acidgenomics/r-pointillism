#' Cell-type markers
#'
#' @name CellTypeMarkers
#' @note Updated 2020-02-20.
#'
#' @inheritParams acidroxygen::params
#' @inheritParams basejump::makeGene2SymbolFromEnsembl
#'
#' @return `CellTypeMarkers`.
#'
#' @examples
#'
#'
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

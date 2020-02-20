#' Cell-type markers
#'
#' @export
#' @note Updated 2020-02-20.
#'
#' @inheritParams acidroxygen::params
#' @inheritParams basejump::makeGene2SymbolFromEnsembl
#'
#' @return `CellTypeMarkers`.
#'
#' @examples
CellTypeMarkers <- function(
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
    class <- "CellTypeMarkers"
    data <- .CellMarkers(
        object = object,
        gene2symbol = gene2symbol,
        class = class
    )
    new(Class = class, data)
}

#' @inherit CellTypeMarkers-class
#' @inheritParams basejump.globals::params
#' @export
CellTypeMarkers <- function(object, gene2symbol) {
    class <- "CellTypeMarkers"
    data <- .cellMarkers(
        object = object,
        gene2symbol = gene2symbol,
        class = class
    )
    new(Class = class, data)
}

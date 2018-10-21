#' @inherit CellTypeMarkers-class
#' @inheritParams general
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

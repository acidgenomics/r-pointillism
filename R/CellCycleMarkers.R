#' @inherit CellCycleMarkers-class
#' @inheritParams general
#' @export
CellCycleMarkers <- function(object, gene2symbol) {
    class <- "CellCycleMarkers"
    data <- .cellMarkers(
        object = object,
        gene2symbol = gene2symbol,
        class = class
    )
    new(Class = class, data)
}

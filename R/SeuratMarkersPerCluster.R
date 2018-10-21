#' @rdname SeuratMarkers
#' @export
SeuratMarkersPerCluster <- function(
    object,
    ranges,
    alpha = 0.05
) {
    data <- .seuratMarkers(
        object = object,
        ranges = ranges,
        alpha = alpha
    )
    out <- split(x = data, f = data[["cluster"]], drop = FALSE)
    names(out) <- paste0("cluster", names(out))
    metadata(out) <- metadata(data)
    new(Class = "MarkersPerCluster", out)
}

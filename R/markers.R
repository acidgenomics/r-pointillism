#' @inherit Markers-class
#' @export
#'
#' @description
#' This generator function is designed to take the original return from a Seurat
#' marker analysis and add corresponding gene annotations.
#'
#' @note For [Seurat::FindAllMarkers()] return, rownames are correctly returned
#'   in the `gene` column.
#'
#' @inheritParams general
#' @param object `data.frame`. Unmodified [Seurat::FindMarkers()] or
#'   [Seurat::FindAllMarkers()] return.
#' @param ranges `GRanges`. Gene annotations. Names must correspond to the
#'   rownames. The function will automatically subset the ranges and arrange
#'   them alphabetically.
#'
#' @examples
#' data(seurat_small)
#' object <- seurat_small
#' ranges <- rowRanges(object)
#'
#' ## `FindMarkers()` return.
#' invisible(capture.output(
#'     markers <- Seurat::FindMarkers(
#'         object = object,
#'         ident.1 = "1",
#'         ident.2 = "0"
#'     )
#' ))
#' x <- MarkersPerCluster(object = markers, ranges = ranges)
#' summary(ident_3_sanitized)
#'
#' ## `FindAllMarkers()` return.
#' invisible(capture.output(
#'     markers <- Seurat::FindAllMarkers(object)
#' ))
#' x <- MarkersPerCluster(object = markers, ranges = ranges)
#' summary(x)
Markers <- function(
    object,
    ranges,
    alpha = 0.05
) {
    data <- .markers(
        object = object,
        ranges = ranges,
        alpha = alpha
    )
    rownames(data) <- data[["name"]]
    data[["name"]] <- NULL
    new(Class = "Markers", data)
}

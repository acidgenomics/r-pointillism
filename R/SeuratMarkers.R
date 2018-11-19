# FIXME This needs to error if the input data.frame contains `cluster` column.
# TODO Consider only using `SeuratMarkers()` as a single generator but
# returning `SeuratMarkers` or `SeuratMarkersPerCluster` automatically.



#' @inherit SeuratMarkers-class
#' @export
#'
#' @description
#' This generator function is designed to take the original return from a Seurat
#' marker analysis and add corresponding gene annotations.
#'
#' @note For [Seurat::FindAllMarkers()] return, rownames are correctly returned
#'   in the `gene` column.
#'
#' @inheritParams basejump::params
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
#'         ident.2 = NULL
#'     )
#' ))
#' x <- SeuratMarkers(object = markers, ranges = ranges)
#' summary(x)
#'
#' ## `FindAllMarkers()` return.
#' invisible(capture.output(
#'     markers <- Seurat::FindAllMarkers(object)
#' ))
#' x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
#' summary(x)
SeuratMarkers <- function(
    object,
    ranges,
    alpha = 0.05
) {
    data <- .seuratMarkers(
        object = object,
        ranges = ranges,
        alpha = alpha
    )
    rownames(data) <- data[["name"]]
    data[["name"]] <- NULL
    new(Class = "SeuratMarkers", data)
}



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
    new(Class = "SeuratMarkersPerCluster", out)
}

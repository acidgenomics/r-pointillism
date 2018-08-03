# nocov start



#' Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return [.Deprecated()].
NULL



#' Defunct Functions
#'
#' @name defunct
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return [.Defunct()].
NULL



# v0.1.0 =======================================================================
#' @rdname deprecated
#' @export
fetchPCAData <- function(object, ...) {
    .Defunct(
        new = "fetchReducedDimData",
        package = "pointillism"
    )
}



#' @rdname deprecated
#' @export
fetchTSNEData <- function(object, ...) {
    .Deprecated("fetchReducedDimData")
    fetchReducedDimData(
        object = object,
        reducedDim = "TSNE",
        ...
    )
}



#' @rdname deprecated
#' @export
fetchTSNEExpressionData <- function(object, ...) {
    .Deprecated("fetchReducedDimExpressionData")
    fetchReducedDimExpressionData(
        object = object,
        reducedDim = "TSNE",
        ...
    )
}



#' @rdname deprecated
#' @export
fetchUMAPData <- function(object, ...) {
    .Deprecated("fetchReducedDimData")
    fetchReducedDimData(
        object = object,
        reducedDim = "UMAP",
        ...
    )
}



#' @rdname deprecated
#' @export
fetchUMAPExpressionData <- function(object, ...) {
    .Deprecated("fetchReducedDimExpressionData")
    fetchReducedDimExpressionData(
        object = object,
        reducedDim = "UMAP",
        ...
    )
}



#' @rdname deprecated
#' @export
plotFeatureTSNE <- function(object, ...) {
    .Deprecated("plotFeature")
    plotFeature(
        object = object,
        reducedDim = "TSNE",
        ...
    )
}



#' @rdname deprecated
#' @export
plotFeatureUMAP <- function(object, ...) {
    .Deprecated("plotFeature")
    plotFeature(
        object = object,
        reducedDim = "UMAP",
        ...
    )
}



#' @rdname defunct
#' @export
plotKnownMarkers <- function(...) {
    .Defunct("plotKnownMarkersDetected")
}



#' @rdname deprecated
#' @export
plotMarkerTSNE <- function(object, ...) {
    .Deprecated("plotMarker")
    plotMarker(
        object = object,
        reducedDim = "TSNE",
        ...
    )
}



#' @rdname deprecated
#' @export
plotMarkerUMAP <- function(object, ...) {
    .Deprecated("plotMarker")
    plotMarker(
        object = object,
        reducedDim = "UMAP",
        ...
    )
}



#' @rdname defunct
#' @export
plotTSNEExpressionData <- function(...) {
    .Defunct("plotMarker")
}



#' @rdname deprecated
#' @export
readCellTypeMarkersFile <- function(...) {
    .Deprecated("readCellTypeMarkers")
    readCellTypeMarkers(...)
}



#' @rdname defunct
#' @export
readMarkers <- function(...) {
    .Defunct("readCellTypeMarkers")
}



#' @rdname defunct
#' @export
readMarkersFile <- function(...) {
    .Defunct("readCellTypeMarkers")
}



#' @rdname deprecated
#' @export
sanitizeMarkers <- function(...) {
    .Deprecated("sanitizeSeuratMarkers")
    sanitizeSeuratMarkers(...)
}



# nocov end

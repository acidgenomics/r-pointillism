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



# v0.2.0 =======================================================================
#' @rdname defunct
#' @export
knownMarkersDetected <- function(...) {
    .Defunct()
}



#' @rdname deprecated
#' @export
plotKnownMarkersDetected <- function(...) {
    .Defunct("plotKnownMarkers")
    plotKnownMarkers(...)
}



#' @rdname deprecated
#' @export
readCellTypeMarkers <- function(...) {
    .Defunct("CellTypeMarkers")
}



#' @rdname deprecated
#' @export
sanitizeSeuratMarkers <- function(...) {
    .Deprecated("SeuratMarkers")
    SeuratMarkers(...)
}



# nocov end

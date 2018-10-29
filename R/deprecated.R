# nocov start



#' Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @inheritParams basejump.globals::params
#'
#' @return [.Deprecated()].
NULL



#' Defunct Functions
#'
#' @name defunct
#' @keywords internal
#'
#' @inheritParams basejump.globals::params
#'
#' @return [.Defunct()].
NULL



# v0.2.0 =======================================================================
#' @rdname defunct
#' @export
knownMarkers <- function(...) {
    .Defunct("KnownMarkers")
}



#' @rdname defunct
#' @export
knownMarkersDetected <- function(...) {
    .Defunct("KnownMarkers")
}



#' @rdname defunct
#' @export
plotKnownMarkersDetected <- function(...) {
    .Defunct("plotKnownMarkers")
}



#' @rdname defunct
#' @export
readCellTypeMarkers <- function(...) {
    .Defunct("CellTypeMarkers")
}



#' @rdname defunct
#' @export
sanitizeMarkers <- function(...) {
    .Defunct("SeuratMarkers")
}



# nocov end

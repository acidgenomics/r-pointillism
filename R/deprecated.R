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
    .Defunct("KnownSeuratMarkers")
}



#' @rdname deprecated
#' @export
readCellTypeMarkers <- function(...) {
    .Deprecated("CellTypeMarkers")
    CellTypeMarkers(...)
}



#' @rdname deprecated
#' @export
sanitizeSeuratMarkers <- function(...) {
    .Deprecated("SeuratMarkers")
    SeuratMarkers(...)
}



# nocov end

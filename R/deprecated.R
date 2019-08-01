## nocov start



#' @name deprecated
#' @inherit basejump::deprecated
#' @keywords internal
NULL

#' @name defunct
#' @inherit basejump::defunct
#' @keywords internal
NULL



## v0.2.0 ======================================================================
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



## v0.2.5 ======================================================================
#' @rdname defunct
#' @export
importCellCycleMarkersFromGoogle <- function(...) {
    .Defunct("CellCycleMarkers")
}


#' @rdname defunct
#' @export
importCellTypeMarkersFromGoogle <- function(...) {
    .Defunct("CellTypeMarkers")
}

#' @rdname defunct
#' @export
runZinbwave <- function(...) {
    .Defunct()
}



## nocov end

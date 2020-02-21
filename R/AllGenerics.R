#' @rdname KnownMarkers
#' @export
setGeneric(
    name = "KnownMarkers",
    def = function(markers, known, ...) {
        standardGeneric("KnownMarkers")
    }
)



#' @rdname SeuratMarkers
#' @export
setGeneric(
    name = "SeuratMarkers",
    def = function(object, ...) {
        standardGeneric("SeuratMarkers")
    }
)



#' @rdname SeuratMarkers
#' @export
setGeneric(
    name = "SeuratMarkersPerCluster",
    def = function(object, ...) {
        standardGeneric("SeuratMarkersPerCluster")
    }
)



#' @rdname runSeurat
#' @export
setGeneric(
    name = "runSeurat",
    def = function(object, ...) {
        standardGeneric("runSeurat")
    }
)



#' @rdname runZinbwave
#' @export
setGeneric(
    name = "runZinbwave",
    def = function(Y, ...) {
        standardGeneric("runZinbwave")
    }
)

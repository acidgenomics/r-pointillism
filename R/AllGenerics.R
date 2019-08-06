#' @rdname KnownMarkers-class
#' @export
setGeneric(
    name = "KnownMarkers",
    def = function(markers, known, ...) {
        standardGeneric("KnownMarkers")
    }
)

#' @rdname SeuratMarkers-class
#' @export
setGeneric(
    name = "SeuratMarkers",
    def = function(object, ...) {
        standardGeneric("SeuratMarkers")
    }
)

#' @rdname SeuratMarkersPerCluster-class
#' @export
setGeneric(
    name = "SeuratMarkersPerCluster",
    def = function(object, ...) {
        standardGeneric("SeuratMarkersPerCluster")
    }
)

#' @rdname normalize
#' @export
setGeneric(
    name = "normalize",
    def = function(object, ...) {
        standardGeneric("normalize")
    }
)

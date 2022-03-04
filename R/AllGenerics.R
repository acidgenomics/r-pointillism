#' @rdname KnownMarkers
#' @name KnownMarkers
#' @importFrom AcidGenerics KnownMarkers
#' @usage KnownMarkers(object, ...)
#' @export
NULL

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

#' @rdname plotPCElbow
#' @name plotPCElbow
#' @importFrom AcidGenerics plotPCElbow
#' @usage plotPCElbow(object, ...)
#' @export
NULL

#' @rdname plotTopMarkers
#' @name plotTopMarkers
#' @importFrom AcidGenerics plotTopMarkers
#' @usage plotTopMarkers(object, markers, ...)
#' @export
NULL

#' @rdname runSeurat
#' @export
setGeneric(
    name = "runSeurat",
    def = function(object, ...) {
        standardGeneric("runSeurat")
    }
)

#' @export
#' @name summary
#' @rdname summary
#' @usage summary(object, ...)
NULL

#' @rdname topMarkers
#' @name topMarkers
#' @importFrom AcidGenerics topMarkers
#' @usage topMarkers(object, ...)
#' @export
NULL

#' @rdname KnownMarkers
#' @name KnownMarkers
#' @importFrom AcidSingleCell KnownMarkers
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



#' @rdname runZinbwave
#' @export
setGeneric(
    name = "runZinbwave",
    def = function(Y, ...) {  # nolint
        standardGeneric("runZinbwave")
    }
)



#' @rdname summary
#' @name summary
#' @importFrom AcidGenerics summary
#' @usage summary(object, ...)
#' @export
NULL



#' @rdname topMarkers
#' @name topMarkers
#' @importFrom AcidGenerics topMarkers
#' @usage topMarkers(object, ...)
#' @export
NULL

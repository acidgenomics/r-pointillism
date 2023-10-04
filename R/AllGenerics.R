#' @rdname KnownMarkers
#' @name KnownMarkers
#' @usage KnownMarkers(markers, known, ...)
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

#' @rdname as.Seurat
#' @name as.Seurat
#' @usage as.Seurat(x, ...)
#' @export
NULL

#' @rdname as.SingleCellExperiment
#' @name as.SingleCellExperiment
#' @usage as.SingleCellExperiment(x, ...)
#' @export
NULL

#' @rdname plotPcElbow
#' @name plotPcElbow
#' @usage plotPcElbow(object, ...)
#' @export
NULL

#' @rdname plotTopMarkers
#' @name plotTopMarkers
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
#' @usage topMarkers(object, ...)
#' @export
NULL

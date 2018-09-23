#' @rdname cellCountsPerCluster
#' @export
setGeneric(
    name = "cellCountsPerCluster",
    def = function(object, ...) {
        standardGeneric("cellCountsPerCluster")
    }
)



#' @rdname cellTypesPerCluster
#' @export
setGeneric(
    name = "cellTypesPerCluster",
    def = function(object, ...) {
        standardGeneric("cellTypesPerCluster")
    }
)



#' @rdname clusterCellCountsPerSample
#' @export
setGeneric(
    name = "clusterCellCountsPerSample",
    def = function(object, ...) {
        standardGeneric("clusterCellCountsPerSample")
    }
)



#' @rdname clusterID
#' @export
setGeneric(
    name = "clusterID",
    def = function(object, ...) {
        standardGeneric("clusterID")
    }
)



#' @rdname diffExp
#' @export
setGeneric(
    name = "diffExp",
    def = function(object, ...) {
        standardGeneric("diffExp")
    }
)



#' @rdname diffExpPerCluster
#' @export
setGeneric(
    name = "diffExpPerCluster",
    def = function(object, ...) {
        standardGeneric("diffExpPerCluster")
    }
)



#' @rdname findMarkers
#' @export
setGeneric(
    name = "findMarkers",
    def = function(object, ...) {
        standardGeneric("findMarkers")
    }
)



#' @rdname knownMarkers
#' @export
setGeneric(
    name = "knownMarkers",
    def = function(all, known, ...) {
        standardGeneric("knownMarkers")
    }
)



#' @rdname plotCellTypesPerCluster
#' @export
setGeneric(
    name = "plotCellTypesPerCluster",
    def = function(object, ...) {
        standardGeneric("plotCellTypesPerCluster")
    }
)



#' @rdname plotGene
#' @export
setGeneric(
    name = "plotDot",
    def = function(object, ...) {
        standardGeneric("plotDot")
    }
)



#' @rdname plotFeature
#' @export
setGeneric(
    name = "plotFeature",
    def = function(object, ...) {
        standardGeneric("plotFeature")
    }
)



#' @rdname plotKnownMarkers
#' @export
setGeneric(
    name = "plotKnownMarkers",
    def = function(object, ...) {
        standardGeneric("plotKnownMarkers")
    }
)



#' @rdname plotMarker
#' @export
setGeneric(
    name = "plotMarker",
    def = function(object, ...) {
        standardGeneric("plotMarker")
    }
)



#' @rdname plotPCElbow
#' @export
setGeneric(
    name = "plotPCElbow",
    def = function(object, ...) {
        standardGeneric("plotPCElbow")
    }
)



#' @rdname plotReducedDim
#' @export
setGeneric(
    name = "plotReducedDim",
    def = function(object, ...) {
        standardGeneric("plotReducedDim")
    }
)



#' @rdname plotTopMarkers
#' @export
setGeneric(
    name = "plotTopMarkers",
    def = function(object, ...) {
        standardGeneric("plotTopMarkers")
    }
)



#' @rdname plotReducedDim
#' @export
setGeneric(
    name = "plotTSNE",
    def = function(object, ...) {
        standardGeneric("plotTSNE")
    }
)



#' @rdname plotReducedDim
#' @export
setGeneric(
    name = "plotUMAP",
    def = function(object, ...) {
        standardGeneric("plotUMAP")
    }
)



#' @rdname plotGene
#' @export
setGeneric(
    name = "plotViolin",
    def = function(object, ...) {
        standardGeneric("plotViolin")
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



#' @rdname topMarkers
#' @export
setGeneric(
    name = "topMarkers",
    def = function(object, ...) {
        standardGeneric("topMarkers")
    }
)

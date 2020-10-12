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



#' @rdname cellCountsPerCluster
#' @name cellCountsPerCluster
#' @importFrom AcidGenerics cellCountsPerCluster
#' @usage cellCountsPerCluster(object, ...)
#' @export
NULL



#' @rdname cellTypesPerCluster
#' @name cellTypesPerCluster
#' @importFrom AcidGenerics cellTypesPerCluster
#' @usage cellTypesPerCluster(object, ...)
#' @export
NULL



#' @rdname clusters
#' @name clusters
#' @importFrom AcidGenerics clusters
#' @usage clusters(object, ...)
#' @export
NULL



#' @rdname cpm
#' @name cpm
#' @importFrom AcidGenerics cpm
#' @usage cpm(object, ...)
#' @export
NULL



#' @rdname diffExp
#' @name diffExp
#' @importFrom AcidGenerics diffExp
#' @usage diffExp(object, ...)
#' @export
NULL



#' @rdname diffExpPerCluster
#' @name diffExpPerCluster
#' @importFrom AcidGenerics diffExpPerCluster
#' @usage diffExpPerCluster(object, ...)
#' @export
NULL



#' @rdname findMarkers
#' @name findMarkers
#' @importFrom AcidGenerics findMarkers
#' @usage findMarkers(object, ...)
#' @export
NULL



#' @rdname normalize
#' @name normalize
#' @importFrom BiocGenerics normalize
#' @usage normalize(object, ...)
#' @export
NULL



#' @rdname plotCellCountsPerCluster
#' @name plotCellCountsPerCluster
#' @importFrom AcidGenerics plotCellCountsPerCluster
#' @usage plotCellCountsPerCluster(object, ...)
#' @export
NULL



#' @rdname plotCellTypesPerCluster
#' @name plotCellTypesPerCluster
#' @importFrom AcidGenerics plotCellTypesPerCluster
#' @usage plotCellTypesPerCluster(object, markers, ...)
#' @export
NULL



#' @rdname plotCounts
#' @name plotCounts
#' @importFrom AcidGenerics plotCounts
#' @usage plotCounts(object, ...)
#' @export
NULL

#' @rdname plotCounts
#' @name plotDots
#' @importFrom AcidGenerics plotDots
#' @usage plotDots(object, ...)
#' @export
NULL

#' @rdname plotCounts
#' @name plotViolin
#' @importFrom AcidGenerics plotViolin
#' @usage plotViolin(object, ...)
#' @export
NULL



#' @rdname plotFeature
#' @name plotFeature
#' @importFrom AcidGenerics plotFeature
#' @usage plotFeature(object, ...)
#' @export
NULL



#' @rdname plotKnownMarkers
#' @name plotKnownMarkers
#' @importFrom AcidGenerics plotKnownMarkers
#' @usage plotKnownMarkers(object, markers, ...)
#' @export
NULL



#' @rdname plotMarker
#' @name plotMarker
#' @importFrom AcidGenerics plotMarker
#' @usage plotMarker(object, ...)
#' @export
NULL



#' @rdname plotPCElbow
#' @name plotPCElbow
#' @importFrom AcidGenerics plotPCElbow
#' @usage plotPCElbow(object, ...)
#' @export
NULL



#' @rdname plotReducedDim
#' @name plotReducedDim
#' @importFrom AcidGenerics plotReducedDim
#' @usage plotReducedDim(object, ...)
#' @export
NULL

#' @rdname plotReducedDim
#' @name plotPCA
#' @importFrom BiocGenerics plotPCA
#' @usage plotPCA(object, ...)
#' @export
NULL

#' @rdname plotReducedDim
#' @name plotTSNE
#' @importFrom AcidGenerics plotTSNE
#' @usage plotTSNE(object, ...)
#' @export
NULL

#' @rdname plotReducedDim
#' @name plotUMAP
#' @importFrom AcidGenerics plotUMAP
#' @usage plotUMAP(object, ...)
#' @export
NULL



#' @rdname plotStackedBarPlot
#' @name plotStackedBarPlot
#' @importFrom AcidGenerics plotStackedBarPlot
#' @usage plotStackedBarPlot(object, ...)
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
#' @importFrom S4Vectors summary
#' @usage summary(object, ...)
#' @export
NULL



#' @rdname topMarkers
#' @name topMarkers
#' @importFrom AcidGenerics topMarkers
#' @usage topMarkers(object, ...)
#' @export
NULL

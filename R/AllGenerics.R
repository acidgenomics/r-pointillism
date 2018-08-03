#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioSingleCell barcodeRanksPerSample
#' @export
setGeneric(
    "barcodeRanksPerSample",
    getGeneric(
        "barcodeRanksPerSample",
        package = "bcbioSingleCell"
    )
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioSingleCell cell2sample
#' @export
setGeneric(
    "cell2sample",
    getGeneric(
        "cell2sample",
        package = "bcbioSingleCell"
    )
)



#' @rdname cellCountsPerCluster
#' @export
setGeneric(
    "cellCountsPerCluster",
    function(object, ...) {
        standardGeneric("cellCountsPerCluster")
    }
)



#' @rdname cellTypesPerCluster
#' @export
setGeneric(
    "cellTypesPerCluster",
    function(object, ...) {
        standardGeneric("cellTypesPerCluster")
    }
)



#' @rdname clusterCellCountsPerSample
#' @export
setGeneric(
    "clusterCellCountsPerSample",
    function(object, ...) {
        standardGeneric("clusterCellCountsPerSample")
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump convertGenesToSymbols
#' @export
setGeneric(
    "convertGenesToSymbols",
    getGeneric(
        "convertGenesToSymbols",
        package = "basejump"
    )
)



#' @rdname diffExp
#' @export
setGeneric(
    "diffExp",
    function(object, ...) {
        standardGeneric("diffExp")
    }
)



#' @rdname fetchData
#' @export
setGeneric(
    "fetchGeneData",
    function(object, ...) {
        standardGeneric("fetchGeneData")
    }
)



#' @rdname fetchData
#' @export
setGeneric(
    "fetchReducedDimData",
    function(object, ...) {
        standardGeneric("fetchReducedDimData")
    }
)



#' @rdname fetchData
#' @export
setGeneric(
    "fetchReducedDimExpressionData",
    function(object, ...) {
        standardGeneric("fetchReducedDimExpressionData")
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump gene2symbol
#' @export
setGeneric(
    "gene2symbol",
    getGeneric(
        "gene2symbol",
        package = "basejump"
    )
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups
#' @export
setGeneric(
    "interestingGroups",
    getGeneric(
        "interestingGroups",
        package = "basejump"
    )
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups<-
#' @export
setGeneric(
    "interestingGroups<-",
    getGeneric(
        "interestingGroups<-",
        package = "basejump"
    )
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioBase metrics
#' @export
setGeneric(
    "metrics",
    getGeneric(
        "metrics",
        package = "bcbioBase"
    )
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioSingleCell metricsPerSample
#' @export
setGeneric(
    "metricsPerSample",
    getGeneric(
        "metricsPerSample",
        package = "bcbioSingleCell"
    )
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioSingleCell plotBarcodeRanks
#' @export
setGeneric(
    "plotBarcodeRanks",
    getGeneric(
        "plotBarcodeRanks",
        package = "bcbioSingleCell"
    )
)



#' @rdname plotCellTypesPerCluster
#' @export
setGeneric(
    "plotCellTypesPerCluster",
    function(object, ...) {
        standardGeneric("plotCellTypesPerCluster")
    }
)



#' @rdname plotFeature
#' @export
setGeneric(
    "plotFeature",
    function(object, ...) {
        standardGeneric("plotFeature")
    }
)



#' @rdname plotMarker
#' @export
setGeneric(
    "plotKnownMarkersDetected",
    function(object, ...) {
        standardGeneric("plotKnownMarkersDetected")
    }
)



#' @rdname plotMarker
#' @export
setGeneric(
    "plotMarker",
    function(object, ...) {
        standardGeneric("plotMarker")
    }
)



#' @rdname plotPCElbow
#' @export
setGeneric(
    "plotPCElbow",
    function(object, ...) {
        standardGeneric("plotPCElbow")
    }
)



#' @rdname plotReducedDim
#' @export
setGeneric(
    "plotReducedDim",
    function(object, ...) {
        standardGeneric("plotReducedDim")
    }
)



#' @rdname plotMarker
#' @export
setGeneric(
    "plotTopMarkers",
    function(object, ...) {
        standardGeneric("plotTopMarkers")
    }
)



#' @rdname plotReducedDim
#' @export
setGeneric(
    "plotTSNE",
    function(object, ...) {
        standardGeneric("plotTSNE")
    }
)



#' @rdname plotReducedDim
#' @export
setGeneric(
    "plotUMAP",
    function(object, ...) {
        standardGeneric("plotUMAP")
    }
)

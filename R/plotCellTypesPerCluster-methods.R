## FIXME This isn't working as expected.

## Error in result[[njob]] <- value :
##   attempt to select less than one element in OneIndex
## Calls: plotCellTypesPerCluster ... .local -> bplapply -> bplapply -> bploop -> bploop.lapply
## Backtrace:
##     █
##  1. ├─bioverbs::plotCellTypesPerCluster(...)
##  2. └─pointillism::plotCellTypesPerCluster(...)
##  3.   └─pointillism:::.local(object, markers, ...)
##  4.     ├─BiocParallel::bplapply(...) R/plotCellTypesPerCluster-methods.R:66:8
##  5.     └─BiocParallel::bplapply(...)
##  6.       ├─BiocParallel::bploop(...)
##  7.       └─BiocParallel:::bploop.lapply(...)



#' @name plotCellTypesPerCluster
#' @inherit bioverbs::plotCellTypesPerCluster
#' @note Updated 2019-09-03.
#'
#' @details
#' Plot the geometric mean of the significant marker genes for every known cell
#' type (per unbiased cluster). Cell types with too few (`min` cutoff) or too
#' many (`max` cutoff) marker genes will be skipped.
#'
#' @inheritParams acidroxygen::params
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @return Show graphical output. Invisibly return `list`.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#' data(seuratKnownMarkers)
#'
#' ## Seurat ====
#' object <- Seurat
#' markers <- seuratKnownMarkers
#'
#' plotCellTypesPerCluster(
#'     object = object,
#'     markers = markers,
#'     reduction = "UMAP"
#' )
NULL



#' @rdname plotCellTypesPerCluster
#' @name plotCellTypesPerCluster
#' @importFrom bioverbs plotCellTypesPerCluster
#' @usage plotCellTypesPerCluster(object, markers, ...)
#' @export
NULL



## Updated 2019-09-03.
`plotCellTypesPerCluster,SingleCellExperiment,KnownMarkers` <-  # nolint
    function(
        object,
        markers,
        min = 1L,
        max = Inf,
        reduction,
        expression,
        headerLevel = 2L,
        ...
    ) {
        ## Passthrough: color, dark.
        validObject(object)
        validObject(markers)
        assert(
            isScalar(reduction),
            isHeaderLevel(headerLevel)
        )
        expression <- match.arg(expression)
        markers <- cellTypesPerCluster(markers, min = min, max = max)
        ## Output Markdown headers per cluster.
        clusters <- unique(as.character(markers[["cluster"]]))
        assert(isNonEmpty(clusters))
        return <- lapply(
            X = clusters,
            FUN = function(cluster) {
                markdownHeader(
                    text = paste("Cluster", cluster),
                    level = headerLevel,
                    tabset = TRUE,
                    asis = TRUE
                )
                keep <- markers[["cluster"]] == cluster
                clusterData <- markers[keep, , drop = FALSE]
                if (!hasRows(clusterData)) {
                    message(sprintf("No markers for cluster %s.", cluster))
                    return(invisible())
                }
                cellTypes <- clusterData[["cellType"]]
                assert(is.factor(cellTypes))
                lapply(
                    X = cellTypes,
                    FUN = function(cellType) {
                        title <- as.character(cellType)
                        markdownHeader(
                            text = title,
                            level = headerLevel + 1L,
                            asis = TRUE
                        )
                        ## Modify the title by adding the cluster number.
                        title <- paste(paste0("Cluster ", cluster, ":"), title)
                        cellData <- clusterData[
                            clusterData[["cellType"]] == cellType,
                            ,
                            drop = FALSE
                            ]
                        assert(identical(nrow(cellData), 1L))
                        genes <- cellData[["name"]]
                        genes <- as.character(genes)
                        genes <- strsplit(genes, split = ", ")[[1L]]
                        if (!hasLength(genes)) {
                            return(invisible())
                        }
                        p <- plotMarker(
                            object = object,
                            genes = genes,
                            reduction = reduction,
                            expression = expression,
                            ...
                        )
                        show(p)
                        invisible(p)
                    }
                )
            }
        )
        invisible(return)
    }

formals(`plotCellTypesPerCluster,SingleCellExperiment,KnownMarkers`)[
    c("reduction", "expression", "BPPARAM")
] <- list(reduction, expression, BPPARAM)



#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    f = "plotCellTypesPerCluster",
    signature = signature(
        object = "SingleCellExperiment",
        markers = "KnownMarkers"
    ),
    definition = `plotCellTypesPerCluster,SingleCellExperiment,KnownMarkers`
)



## Updated 2019-08-07.
`plotCellTypesPerCluster,Seurat,KnownMarkers` <-  # nolint
    `plotCellTypesPerCluster,SingleCellExperiment,KnownMarkers`



#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    f = "plotCellTypesPerCluster",
    signature = signature(
        object = "Seurat",
        markers = "KnownMarkers"
    ),
    definition = `plotCellTypesPerCluster,Seurat,KnownMarkers`
)

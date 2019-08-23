#' @name plotCellTypesPerCluster
#' @inherit bioverbs::plotCellTypesPerCluster
#' @note Updated 2019-08-23.
#'
#' @details
#' Plot the geometric mean of the significant marker genes for every known cell
#' type (per unbiased cluster). Cell types with too few (`min` cutoff) or too
#' many (`max` cutoff) marker genes will be skipped.
#'
#' @inheritParams acidroxygen::params
#' @inheritParams BiocParallel::bplapply
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



## Updated 2019-08-23.
`plotCellTypesPerCluster,SingleCellExperiment,KnownMarkers` <-  # nolint
    function(
        object,
        markers,
        min = 1L,
        max = Inf,
        reduction,
        expression,
        headerLevel = 2L,
        BPPARAM = BiocParallel::SerialParam(),  # nolint
        ...
    ) {
        ## Passthrough: color, dark.
        validObject(object)
        validObject(markers)
        assert(isScalar(reduction))
        expression <- match.arg(expression)
        assert(isHeaderLevel(headerLevel))
        markers <- cellTypesPerCluster(
            object = markers,
            min = min,
            max = max
        )
        assert(
            is(markers, "grouped_df"),
            hasRows(markers)
        )
        ## Output Markdown headers per cluster.
        clusters <- markers[["cluster"]] %>%
            as.character() %>%
            unique()
        assert(isNonEmpty(clusters))
        return <- bplapply(
            X = clusters,
            FUN = function(cluster) {
                markdownHeader(
                    text = paste("Cluster", cluster),
                    level = headerLevel,
                    tabset = TRUE,
                    asis = TRUE
                )
                clusterData <- filter(markers, !!sym("cluster") == !!cluster)
                if (nrow(clusterData) == 0L) {
                    message(sprintf("No markers for cluster %s.", cluster))
                    return(invisible())
                }
                assert(hasRows(clusterData))
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
                        cellData <-
                            filter(clusterData, !!sym("cellType") == !!cellType)
                        assert(nrow(cellData) == 1L)
                        genes <- cellData %>%
                            pull("name") %>%
                            as.character() %>%
                            strsplit(", ") %>%
                            .[[1L]]
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
            },
            BPPARAM = BPPARAM
        )
        invisible(return)
    }

formals(`plotCellTypesPerCluster,SingleCellExperiment,KnownMarkers`)[
    c("reduction", "expression")
] <- list(reduction, expression)



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

#' @name plotCellTypesPerCluster
#' @include globals.R
#' @inherit bioverbs::plotCellTypesPerCluster
#' @inheritParams basejump::params
#'
#' @details
#' Plot the geometric mean of the significant marker genes for every known cell
#' type (per unbiased cluster). Cell types with too few (`min` cutoff) or too
#' many (`max` cutoff) marker genes will be skipped.
#'
#' @param markers `KnownMarkers`.
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @return Show graphical output. Invisibly return `list`.
#'
#' @examples
#' data(seurat_small, known_markers_small)
#' plotCellTypesPerCluster(
#'     object = seurat_small,
#'     markers = known_markers_small
#' )
NULL



#' @importFrom bioverbs plotCellTypesPerCluster
#' @aliases NULL
#' @export
bioverbs::plotCellTypesPerCluster



plotCellTypesPerCluster.SingleCellExperiment <-  # nolint
    function(
        object,
        markers,
        min = 1L,
        max = Inf,
        reducedDim,
        expression,
        headerLevel = 2L,
        progress = FALSE,
        ...
    ) {
        # Passthrough: color, dark.
        validObject(object)
        validObject(markers)
        assert(isScalar(reducedDim))
        expression <- match.arg(expression)
        assert(
            isHeaderLevel(headerLevel),
            isFlag(progress)
        )
        if (isTRUE(progress)) {
            applyFun <- pblapply
        } else {
            applyFun <- lapply
        }

        markers <- cellTypesPerCluster(
            object = markers,
            min = min,
            max = max
        )
        assert(
            is(markers, "grouped_df"),
            hasRows(markers)
        )

        # Output Markdown headers per cluster.
        clusters <- markers[["cluster"]] %>%
            as.character() %>%
            unique()
        assert(isNonEmpty(clusters))

        return <- applyFun(clusters, function(cluster) {
            markdownHeader(
                text = paste("Cluster", cluster),
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )
            clusterData <- filter(markers, !!sym("cluster") == !!cluster)
            if (nrow(clusterData) == 0L) {
                message(paste0("No markers for cluster ", cluster, "."))
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
                    # Modify the title by adding the cluster number.
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
                        reducedDim = reducedDim,
                        expression = expression,
                        ...
                    )
                    show(p)
                    invisible(p)
                }
            )
        })

        invisible(return)
    }
formals(plotCellTypesPerCluster.SingleCellExperiment)[["reducedDim"]] <-
    reducedDim
formals(plotCellTypesPerCluster.SingleCellExperiment)[["expression"]] <-
    expression



#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    f = "plotCellTypesPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = plotCellTypesPerCluster.SingleCellExperiment
)



#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    f = "plotCellTypesPerCluster",
    signature = signature("seurat"),
    definition = plotCellTypesPerCluster.SingleCellExperiment
)

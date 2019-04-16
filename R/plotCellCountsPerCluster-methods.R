#' @name plotCellCountsPerCluster
#' @inherit bioverbs::plotCellCountsPerCluster
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @return Show graphical output. Invisibly return `ggplot`.
#'
#' @examples
#' data(seurat)
#' plotCellCountsPerCluster(seurat)
NULL



#' @rdname plotCellCountsPerCluster
#' @name plotCellCountsPerCluster
#' @importFrom bioverbs plotCellCountsPerCluster
#' @export
NULL



plotCellCountsPerCluster.SingleCellExperiment <-  # nolint
    function(
        object,
        interestingGroups = NULL
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

        markers <- CellCountsPerCluster(
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
            CellCounts <- clusterData[["cellType"]]
            assert(is.factor(CellCounts))
            lapply(
                X = CellCounts,
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



#' @rdname plotCellCountsPerCluster
#' @export
setMethod(
    f = "plotCellCountsPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = plotCellCountsPerCluster.SingleCellExperiment
)



plotCellCountsPerCluster.Seurat <-  # nolint
    plotCellCountsPerCluster.SingleCellExperiment



#' @rdname plotCellCountsPerCluster
#' @export
setMethod(
    f = "plotCellCountsPerCluster",
    signature = signature("Seurat"),
    definition = plotCellCountsPerCluster.Seurat
)

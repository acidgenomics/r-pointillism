#' Plot Cell Types per Cluster
#'
#' Plot the geometric mean of the significant marker genes for every known cell
#' type (per unbiased cluster). Cell types with too few (`min` cutoff) or too
#' many (`max` cutoff) marker genes will be skipped.
#'
#' @name plotCellTypesPerCluster
#' @include globals.R
#'
#' @inheritParams basejump::params
#' @param markers `KnownMarkers`.
#' @param ... Passthrough arguments to `plotMarker()`.
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
        assert_is_scalar(reducedDim)
        expression <- match.arg(expression)
        assertIsHeaderLevel(headerLevel)
        assert_is_a_bool(progress)
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
        assert_is_all_of(markers, "grouped_df")
        assert_has_rows(markers)

        # Output Markdown headers per cluster.
        clusters <- markers[["cluster"]] %>%
            as.character() %>%
            unique()
        assert_is_non_empty(clusters)

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
            assert_has_rows(clusterData)
            cellTypes <- clusterData[["cellType"]]
            assert_is_factor(cellTypes)
            lapply(
                X = cellTypes,
                FUN = function(cellType) {
                    cellData <- clusterData %>%
                        filter(!!sym("cellType") == !!cellType)
                    assert_that(nrow(cellData) == 1L)
                    genes <- cellData %>%
                        pull("name") %>%
                        as.character() %>%
                        strsplit(", ") %>%
                        .[[1L]]
                    title <- as.character(cellType)
                    markdownHeader(
                        text = title,
                        level = headerLevel + 1L,
                        asis = TRUE
                    )
                    # Modify the title by adding the cluster number.
                    title <- paste(paste0("Cluster ", cluster, ":"), title)
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
    definition = getMethod(
        f = "plotCellTypesPerCluster",
        signature = signature("SingleCellExperiment")
    )
)

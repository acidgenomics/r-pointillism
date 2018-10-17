# FIXME Fix the formals.
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function 'mapGenesToSymbols' for signature '"seurat"'
# FIXME Need to update to use SeuratMarkers class.
# FIXME Add progress option, switching pbapply.
# FIXME Set reducedDim to "TSNE" by default.
# FIXME Ensure "PC1" in axis label instead of "pc1"



#' Plot Cell Types per Cluster
#'
#' Plot the geometric mean of the significant marker genes for every known cell
#' type (per unbiased cluster). Cell types with too few (`min` cutoff) or too
#' many (`max` cutoff) marker genes will be skipped.
#'
#' @name plotCellTypesPerCluster
#' @family Plot Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param markers `KnownSeuratMarkers`.
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



# FIXME Set these formals automatically.
.plotCellTypesPerCluster.SCE <-  # nolint
    function(
        object,
        markers,
        min = 1L,
        max = Inf,
        reducedDim = 1L,
        expression = c("mean", "median", "sum"),
        headerLevel = 2L
    ) {
        # Passthrough: color, dark.
        validObject(object)
        validObject(markers)
        assert_is_scalar(reducedDim)
        expression <- match.arg(expression)
        assertIsHeaderLevel(headerLevel)

        markers <- cellTypesPerCluster(object = markers, min = min, max = max)
        assert_is_all_of(markers, "grouped_df")
        assert_has_rows(markers)

        # Output Markdown headers per cluster
        clusters <- levels(markers[["cluster"]])
        assert_is_non_empty(clusters)

        return <- pblapply(clusters, function(cluster) {
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
                    stopifnot(nrow(cellData) == 1L)
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
                        expression = expression
                    )
                    show(p)
                    invisible(p)
                }
            )
        })

        invisible(return)
    }



#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    f = "plotCellTypesPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = .plotCellTypesPerCluster.SCE
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

# FIXME
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function 'mapGenesToSymbols' for signature '"seurat"'



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
#' @param markers `grouped_df`. Cell types per cluster data. This must be the
#'   return from [cellTypesPerCluster()].
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @return Show graphical output. Invisibly return `list`.
#'
#' @examples
#' markers <- cellTypesPerCluster(known_markers_detected_small)
#' # Subset for speed.
#' markers <- head(markers, n = 2L)
#' glimpse(markers)
#'
#' # Let's plot the first row, as an example.
#' plotCellTypesPerCluster(
#'     object = seurat_small,
#'     markers = markers
#' )
NULL



#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    "plotCellTypesPerCluster",
    signature("SingleCellExperiment"),
    function(
        object,
        markers,
        reducedDim = c("TSNE", "UMAP"),
        expression = c("mean", "median", "sum"),
        headerLevel = 2L,
        ...
    ) {
        # Legacy arguments -----------------------------------------------------
        call <- match.call()
        if ("cellTypesPerCluster" %in% names(call)) {
            stop("Use `markers` instead of `cellTypesPerCluster`")
        }

        # Assert checks --------------------------------------------------------
        # Passthrough: color, dark
        validObject(object)
        stopifnot(is(markers, "grouped_df"))
        assert_has_rows(markers)
        assert_are_identical(
            x = group_vars(markers),
            y = "cluster"
        )
        assert_is_subset(c("cellType", "rowname"), colnames(markers))
        reducedDim <- match.arg(reducedDim)
        expression <- match.arg(expression)
        assertIsAHeaderLevel(headerLevel)

        markers <- markers %>%
            ungroup() %>%
            mutate_if(is.factor, droplevels)

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
            subset <- markers %>%
                .[.[["cluster"]] == cluster, , drop = FALSE]
            assert_has_rows(subset)
            lapply(seq_len(nrow(subset)), function(x) {
                cellType <- subset[x, , drop = FALSE]
                genes <- pull(cellType, "rowname") %>%
                    as.character() %>%
                    strsplit(", ") %>%
                    .[[1L]]
                title <- as.character(pull(cellType, "cellType"))
                markdownHeader(
                    text = title,
                    level = headerLevel + 1L,
                    asis = TRUE
                )
                # Modify the title by adding the cluster number (for the plot)
                title <- paste(paste0("Cluster ", cluster, ":"), title)
                p <- plotMarker(
                    object = object,
                    genes = genes,
                    reducedDim = reducedDim,
                    expression,
                    ...
                )
                show(p)
                invisible(p)
            })
        })

        invisible(return)
    }
)



#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    "plotCellTypesPerCluster",
    signature("seurat"),
    getMethod("plotCellTypesPerCluster", "SingleCellExperiment")
)

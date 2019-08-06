## FIXME Need to add cell_data_set support here.



#' @name plotTopMarkers
#' @include globals.R
#' @inherit bioverbs::plotTopMarkers
#' @note Updated 2019-07-31.
#'
#' @details
#' The number of markers to plot is determined by the output of the
#' [topMarkers()] function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.

#' @inheritParams acidroxygen::params
#' @inheritParams params
#' @inheritParams topMarkers
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @examples
#' data(Seurat, package = "acidtest")
#' data(seuratAllMarkers)
#'
#' ## Seurat, SeuratMarkersPerCluster ====
#' object <- Seurat
#' markers <- seuratAllMarkers
#' plotTopMarkers(object = object, markers = markers)
NULL



#' @rdname plotTopMarkers
#' @name plotTopMarkers
#' @importFrom bioverbs plotTopMarkers
#' @usage plotTopMarkers(object, ...)
#' @export
NULL



## Updated 2019-08-02.
`plotTopMarkers,Seurat,SeuratMarkersPerCluster` <-  # nolint
    function(
        object,
        markers,
        n = 1L,
        direction,
        reduction,
        headerLevel = 2L,
        progress = FALSE,
        ...
    ) {
        ## Passthrough: n, direction, coding
        validObject(object)
        validObject(markers)
        markers <- topMarkers(
            object = markers,
            n = n,
            direction = direction
        )
        assert(
            isScalar(reduction),
            isHeaderLevel(headerLevel),
            isFlag(progress)
        )
        if (isTRUE(progress)) {
            applyFun <- pblapply
        } else {
            applyFun <- lapply
        }

        assert(isSubset("cluster", colnames(markers)))
        clusters <- levels(markers[["cluster"]])


        list <- applyFun(clusters, function(cluster) {
            genes <- markers %>%
                ungroup() %>%
                filter(cluster == !!cluster) %>%
                pull("name")
            if (!length(genes)) {
                message(paste0("No genes for cluster ", cluster, "."))
                return(invisible())
            }
            if (length(genes) > 10L) {
                warning("Maximum of 10 genes per cluster is recommended")
            }

            markdownHeader(
                text = paste("Cluster", cluster),
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )

            lapply(genes, function(gene) {
                p <- plotMarker(
                    object = object,
                    genes = gene,
                    reduction = reduction,
                    ...
                )
                show(p)
                invisible(p)
            })
        })

        invisible(list)
    }

formals(`plotTopMarkers,Seurat,SeuratMarkersPerCluster`)[
    c("direction", "reduction")
] <- list(direction, reduction)



#' @rdname plotTopMarkers
#' @export
setMethod(
    f = "plotTopMarkers",
    signature = signature(
        object = "Seurat",
        markers = "SeuratMarkersPerCluster"
    ),
    definition = `plotTopMarkers,Seurat,SeuratMarkersPerCluster`
)

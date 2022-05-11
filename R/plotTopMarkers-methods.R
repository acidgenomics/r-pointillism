#' @name plotTopMarkers
#' @inherit AcidGenerics::plotTopMarkers
#' @note Updated 2022-05-11.
#'
#' @details
#' The number of markers to plot is determined by the output of the
#' `topMarkers()` function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.
#'
#' @inheritParams topMarkers
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough arguments to `plotMarker()`.
#'
#' @examples
#' data(seurat, smpc)
#'
#' ## Seurat, SeuratMarkersPerCluster ====
#' object <- seurat
#' markers <- smpc
#' plotTopMarkers(
#'     object = object,
#'     markers = markers,
#'     reduction = "UMAP"
#' )
NULL



## Updated 2021-03-03.
`plotTopMarkers,Seurat,SeuratMarkersPerCluster` <-  # nolint
    function(
        object,
        markers,
        direction,
        reduction,
        n = 1L,
        headerLevel = 2L,
        ...
    ) {
        ## Passthrough: n, direction, coding
        validObject(object)
        validObject(markers)
        markers <- topMarkers(
            object = markers,
            direction = direction,
            n = n
        )
        assert(
            is(markers, "DataFrame"),
            isScalar(reduction),
            isHeaderLevel(headerLevel)
        )
        assert(is.factor(markers[["cluster"]]))
        clusters <- levels(markers[["cluster"]])
        list <- lapply(
            X = clusters,
            FUN = function(cluster) {
                genes <- markers[
                    markers[["cluster"]] == cluster,
                    "name",
                    drop = TRUE
                    ]
                genes <- as.character(genes)
                if (!hasLength(genes)) {
                    alertWarning(sprintf(
                        "No genes for cluster %s.", cluster
                    ))
                    return(invisible(NULL))
                } else if (length(genes) > 10L) {
                    alertWarning(
                        "Maximum of 10 genes per cluster is recommended."
                    )
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
            }
        )
        invisible(list)
    }

args <- c("direction", "reduction", "BPPARAM")
formals(`plotTopMarkers,Seurat,SeuratMarkersPerCluster`)[args] <-
    .formalsList[args]
rm(args)



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

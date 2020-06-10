#' @name plotStackedBarPlot
#' @inherit acidgenerics::plotStackedBarPlot
#' @note Updated 2020-06-10.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#' data(seurat_known_markers)
#'
#' ## Seurat ====
#' object <- Seurat
#' markers <- seurat_known_markers
#'
#' plotCellTypesPerCluster(
#'     object = object,
#'     markers = markers,
#'     reduction = "UMAP"
#' )
NULL



#' @rdname plotStackedBarPlot
#' @name plotStackedBarPlot
#' @importFrom acidgenerics plotStackedBarPlot
#' @usage plotStackedBarPlot(object, ...)
#' @export
NULL



`plotStackedBarPlot,SingleCellExperiment` <-  # nolint
    function(object) {
        print("FIXME")
    }



#' @rdname plotStackedBarPlot
#' @export
setMethod(
    f = "plotCellTypesPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = `plotStackedBarPlot,SingleCellExperiment`
)



`plotStackedBarPlot,Seurat` <-  # nolint
    `plotStackedBarPlot,SingleCellExperiment`



#' @rdname plotStackedBarPlot
#' @export
setMethod(
    f = "plotCellTypesPerCluster",
    signature = signature("Seurat"),
    definition = `plotStackedBarPlot,Seurat`
)

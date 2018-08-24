#' Plot Gene
#'
#' Visualize genes on a dot or violin plot.
#'
#' @name plotGene
#' @family Plot Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump plotGene
#' @export
#'
#' @inheritParams general
#' @param geom `string`. Plot type. Uses [match.arg()] to pick the type.
#'   Currently supports "`dot`" and "`violin`".
#'
#' @seealso
#' - [plotDot()], [Seurat::DotPlot()].
#' - [plotViolin()], [Seurat::VlnPlot()].
#' - [Seurat::RidgePlot()].
#'
#' @return `ggplot`.
#'
#' @examples
#' object <- sce_small
#' genes <- head(rownames(object))
#' glimpse(genes)
#'
#' # Dot
#' plotGene(object, genes = genes, geom = "dot")
#'
#' # Violin
#' plotGene(object, genes = genes, geom = "violin")
NULL



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        geom = c("dot", "violin"),
        legend = getOption("pointillism.legend", TRUE)
    ) {
        geom <- match.arg(geom)
        if (geom == "dot") {
            fun <- plotDot
        } else if (geom == "violin") {
            fun <- plotViolin
        }
        fun(
            object = object,
            genes = genes,
            legend = legend
        )
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "plotGene",
    signature("seurat"),
    getMethod("plotGene", "SingleCellExperiment")
)

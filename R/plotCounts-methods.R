#' @name plotCounts
#' @aliases plotDots plotViolin
#' @inherit acidgenerics::plotCounts
#'
#' @note Dot geom currently only supports logcounts.
#' @note Updated 2020-02-21.
#'
#' @description Visualize genes on a dot or violin plot.
#'
#' @inheritParams acidroxygen::params
#' @inheritParams ggplot2::geom_violin
#' @param colMin `numeric(1)`.
#'   Minimum scaled average expression threshold. Everything smaller will be
#'   set to this.
#' @param colMax `numeric(1)`.
#'   Maximum scaled average expression threshold. Everything larger will be set
#'   to this.
#' @param dotMin `numeric(1)`.
#'   The fraction of cells at which to draw the smallest dot. All cell groups
#'   with less than this expressing the given gene will have no dot drawn.
#' @param dotScale `numeric(1)`.
#'   Scale the size of the points, similar to `cex`.
#' @param geom `character(1)`.
#'   Plot type. Uses [`match.arg()`][base::match.arg] to pick the type.
#'   Currently supports `"dot"` and `"violin"`.
#' @param ... Additional arguments.
#'
#' @seealso
#' - `Seurat::DotPlot()`.
#' - `Seurat::VlnPlot()`.
#' - `Seurat::RidgePlot()`.
#' - `monocle3::plot_genes_violin()`.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#'
#' ## Plotting with either gene IDs or gene names (symbols) works.
#' genes <- head(rownames(object), n = 4L)
#' print(genes)
#'
#' ## Per sample mode enabled.
#' plotCounts(object, genes = genes, perSample = TRUE)
#'
#' ## Per sample mode disabled.
#' plotCounts(object, genes = genes, perSample = FALSE)
NULL



#' @rdname plotCounts
#' @name plotCounts
#' @importFrom acidgenerics plotCounts
#' @usage plotCounts(object, ...)
#' @export
NULL



## Updated 2020-02-21.
`plotCounts,SingleCellExperiment` <-  # nolint
    function(
        object,
        genes,
        assay = c("logcounts", "normcounts"),
        geom = c("violin", "dot"),
        perSample = TRUE,
        legend,
        title = NULL
    ) {
        validObject(object)
        assay <- match.arg(assay)
        geom <- match.arg(geom)
        args <- as.list(sys.call(which = -1L))[-1L]
        args[["geom"]] <- NULL
        if (geom == "dot") {
            assert(identical(assay, "logcounts"))
            args[["assay"]] <- NULL
            what <- plotDots
        } else if (geom == "violin") {
            what <- plotViolin
        }
        do.call(what = what, args = args)
    }

args <- "legend"
formals(`plotCounts,SingleCellExperiment`)[args] <- .formalsList[args]
rm(args)



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotCounts",
    signature = signature("SingleCellExperiment"),
    definition = `plotCounts,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotCounts,Seurat` <-  # nolint
    `plotCounts,SingleCellExperiment`



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotCounts",
    signature = signature("Seurat"),
    definition = `plotCounts,SingleCellExperiment`
)

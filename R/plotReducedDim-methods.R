#' @name plotReducedDim
#' @aliases plotPCA plotTSNE plotUMAP
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#' @inherit bioverbs::plotReducedDim
#' @note Updated 2019-08-02.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @details
#' Colors using `ident` column defined in
#' [`colData()`][SummarizedExperiment::colData] by default.
#'
#' @section Reduction types:
#'
#' - PCA: **P**rincipal **C**omponent **A**nalysis.
#' - t-SNE: **t**-distributed **S**tochastic **N**eighbor **E**mbedding.
#' - UMAP: **U**niform **M**anifold **A**pproximation and **P**rojection.
#'
#' @section UMAP calculation:
#'
#' [UMAP][] calculation in R requires the [Python][] module `umap-learn`.
#' The UMAP module can be loaded in R using [reticulate][].
#'
#' [Python]: https://www.python.org
#' [UMAP]: https://github.com/lmcinnes/umap
#' [reticulate]: https://rstudio.github.io/reticulate/
#'
#' @seealso
#' - `Seurat::DimPlot()`.
#' - `monocle3::plot_cells()`.
#' - [Seurat Mouse Cell Atlas vignette](https://satijalab.org/seurat/mca.html).
#'
#' @examples
#' data(
#'     Seurat,
#'     cell_data_set,
#'     package = "acidtest"
#' )
#'
#' ## Seurat ====
#' object <- Seurat
#' plotReducedDim(object, reduction = "UMAP")
#'
#' ## cell_data_set ====
#' object <- cell_data_set
#' plotReducedDim(object, reduction = "UMAP")
NULL



#' @rdname plotReducedDim
#' @name plotReducedDim
#' @importFrom bioverbs plotReducedDim
#' @usage plotReducedDim(object, ...)
#' @export
NULL

#' @rdname plotReducedDim
#' @name plotPCA
#' @importFrom BiocGenerics plotPCA
#' @usage plotPCA(object, ...)
#' @export
NULL

#' @rdname plotReducedDim
#' @name plotTSNE
#' @importFrom bioverbs plotTSNE
#' @usage plotTSNE(object, ...)
#' @export
NULL

#' @rdname plotReducedDim
#' @name plotUMAP
#' @importFrom bioverbs plotUMAP
#' @usage plotUMAP(object, ...)
#' @export
NULL



## Constructors ================================================================
## Updated 2019-08-02.
`plotReducedDim,SingleCellExperiment` <-  # nolint
    function(
        object,
        reduction,
        dims,
        interestingGroups = NULL,
        color,
        pointSize,
        pointAlpha,
        pointsAsNumbers,
        label,
        labelSize,
        dark,
        legend,
        title = NULL
    ) {
        validObject(object)
        if (!is(object, "SingleCellExperiment")) {
            object <- as(object, "SingleCellExperiment")
        }
        assert(
            .hasClusters(object),
            ## Allow pass in of positional scalar, for looping.
            isScalar(reduction),
            hasLength(dims, n = 2L),
            all(isIntegerish(dims)),
            isCharacter(interestingGroups, nullOK = TRUE),
            isGGScale(color, scale = "discrete", aes = "colour", nullOK = TRUE),
            isNumber(pointSize),
            isNumber(pointAlpha),
            isFlag(pointsAsNumbers),
            isFlag(label),
            isNumber(labelSize),
            isFlag(dark),
            isFlag(legend),
            isString(title, nullOK = TRUE)
        )
        if (is.null(interestingGroups)) {
            interestingGroups <- "ident"
        }

        data <- .fetchReductionData(
            object = object,
            reduction = reduction,
            dims = dims
        )
        assert(
            is(data, "DataFrame"),
            isSubset(
                x = c("x", "y", "centerX", "centerY", "interestingGroups"),
                y = colnames(data)
            )
        )

        ## Check if interesting groups input is supported.
        supported <- bapply(data, is.factor)
        supported <- names(supported)[supported]
        blacklist <- c("interestingGroups", "orig.ident", "sampleID")
        supported <- setdiff(supported, blacklist)
        if (!isSubset(interestingGroups, supported)) {
            setdiff <- setdiff(interestingGroups, supported)
            stop(sprintf(
                fmt = paste0(
                    "%s ",
                    ngettext(
                        n = length(setdiff),
                        msg1 = "interesting group",
                        msg2 = "interesting groups"
                    ),
                    " not defined: %s\n",
                    "Available:\n%s"
                ),
                length(setdiff),
                toString(setdiff, width = 200L),
                printString(supported)
            ))
        }

        if (isString(interestingGroups)) {
            data[["interestingGroups"]] <- data[[interestingGroups]]
        } else {
            data[["interestingGroups"]] <- apply(
                X = data[ , interestingGroups],
                MARGIN = 1L,
                FUN = paste ,
                collapse = ":"
            )
        }

        ## Turn off labeling if there's only 1 cluster.
        if (hasLength(levels(data[["ident"]]), n = 1L)) {
            label = FALSE
        }

        ## Set the x- and y-axis labels (e.g. t_SNE1, t_SNE2). We're setting
        ## this up internally as the first two columns in the data frame.
        axes <- colnames(data)[seq_len(2L)]
        p <- ggplot(
            data = as_tibble(data),
            mapping = aes(
                x = !!sym("x"),
                y = !!sym("y"),
                color = !!sym("interestingGroups")
            )
        ) +
            labs(
                x = axes[[1L]],
                y = axes[[2L]],
                color = paste(interestingGroups, collapse = ":\n")
            )

        if (isTRUE(pointsAsNumbers)) {
            ## Increase the size, if necessary.
            if (pointSize < 4L) {
                message("Increase pointSize to 4.")
                pointSize <- 4L
            }
            p <- p +
                geom_text(
                    mapping = aes(
                        x = !!sym("x"),
                        y = !!sym("y"),
                        label = !!sym("ident"),
                        color = !!sym("interestingGroups")
                    ),
                    alpha = pointAlpha,
                    size = pointSize,
                    show.legend = legend
                )
        } else {
            p <- p +
                geom_point(
                    alpha = pointAlpha,
                    size = pointSize,
                    show.legend = legend
                )
        }

        if (isTRUE(label)) {
            if (isTRUE(dark)) {
                labelColor <- "white"
            } else {
                labelColor <- "black"
            }
            p <- p +
                geom_text(
                    mapping = aes(
                        x = !!sym("centerX"),
                        y = !!sym("centerY"),
                        label = !!sym("ident")
                    ),
                    color = labelColor,
                    size = labelSize,
                    fontface = "bold"
                )
        }

        ## Dark mode.
        if (isTRUE(dark)) {
            p <- p + acid_theme_dark()
        }

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        ## Improve the axis breaks.
        p <- p +
            scale_x_continuous(breaks = pretty_breaks(n = 4L)) +
            scale_y_continuous(breaks = pretty_breaks(n = 4L))

        p
    }

formals(`plotReducedDim,SingleCellExperiment`)[c(
    "color",
    "dark",
    "dims",
    "label",
    "labelSize",
    "legend",
    "pointAlpha",
    "pointSize",
    "pointsAsNumbers",
    "reduction"
)] <- list(
    color = discreteColor,
    dark = dark,
    dims = dims,
    label = label,
    labelSize = labelSize,
    legend = legend,
    pointAlpha = pointAlpha,
    pointSize = pointSize,
    pointsAsNumbers = pointsAsNumbers,
    reduction = reduction
)



## Updated 2019-07-31.
`plotPCA,SingleCellExperiment` <-  # nolint
    function() {
        do.call(
            what = plotReducedDim,
            args = matchArgsToDoCall(
                args = list(
                    object = object,
                    reduction = "PCA"
                )
            )
        )

    }



## Updated 2019-07-31.
`plotTSNE,SingleCellExperiment` <-  # nolint
    function() {
        do.call(
            what = plotReducedDim,
            args = matchArgsToDoCall(
                args = list(
                    object = object,
                    reduction = "TSNE"
                )
            )
        )
    }



## Updated 2019-07-31.
`plotUMAP,SingleCellExperiment` <-  # nolint
    function() {
        do.call(
            what = plotReducedDim,
            args = matchArgsToDoCall(
                args = list(
                    object = object,
                    reduction = "UMAP"
                )
            )
        )
    }



## Formals =====================================================================
## Set the formals for the convenience functions.
f <- formals(`plotReducedDim,SingleCellExperiment`)
f <- f[setdiff(names(f), "reduction")]
formals(`plotPCA,SingleCellExperiment`) <- f
formals(`plotTSNE,SingleCellExperiment`) <- f
formals(`plotUMAP,SingleCellExperiment`) <- f
rm(f)



## Methods =====================================================================
#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("SingleCellExperiment"),
    definition = `plotReducedDim,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotReducedDim,Seurat` <-  # nolint
    `plotReducedDim,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("Seurat"),
    definition = `plotReducedDim,Seurat`
)



## Updated 2019-08-02.
`plotReducedDim,cell_data_set` <-  # nolint
    `plotReducedDim,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("cell_data_set"),
    definition = `plotReducedDim,cell_data_set`
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("SingleCellExperiment"),
    definition = `plotTSNE,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotTSNE,Seurat` <-  # nolint
    `plotTSNE,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("Seurat"),
    definition = `plotTSNE,Seurat`
)



## Updated 2019-08-02.
`plotTSNE,cell_data_set` <-  # nolint
    `plotTSNE,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("cell_data_set"),
    definition = `plotTSNE,cell_data_set`
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature("SingleCellExperiment"),
    definition = `plotUMAP,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotUMAP,Seurat` <-  # nolint
    `plotUMAP,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature("Seurat"),
    definition = `plotUMAP,Seurat`
)



## Updated 2019-08-02.
`plotUMAP,cell_data_set` <-  # nolint
    `plotUMAP,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature("cell_data_set"),
    definition = `plotUMAP,cell_data_set`
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("SingleCellExperiment"),
    definition = `plotPCA,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotPCA,Seurat` <-  # nolint
    `plotPCA,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("Seurat"),
    definition = `plotPCA,Seurat`
)



## Updated 2019-08-02.
`plotPCA,cell_data_set` <-  # nolint
    `plotPCA,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("cell_data_set"),
    definition = `plotPCA,cell_data_set`
)

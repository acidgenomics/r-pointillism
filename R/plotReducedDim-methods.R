#' @name plotReducedDim
#' @aliases plotPCA plotTSNE plotUMAP
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit AcidGenerics::plotReducedDim
#' @note Updated 2021-03-03.
#'
#' @inheritParams AcidRoxygen::params
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
#' data(Seurat, package = "AcidTest")
#'
#' ## Seurat ====
#' object <- Seurat
#' plotReducedDim(object, reduction = "UMAP")
NULL



## Updated 2021-03-03.
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
        labels = list(
            title = NULL,
            subtitle = NULL
        )
    ) {
        validObject(object)
        assert(
            .hasClusters(object),
            ## Allow pass in of positional scalar, for looping.
            isScalar(reduction),
            hasLength(dims, n = 2L),
            all(isIntegerish(dims)),
            isCharacter(interestingGroups, nullOK = TRUE),
            isGGScale(color, scale = "discrete", aes = "color", nullOK = TRUE),
            isNumber(pointSize),
            isNumber(pointAlpha),
            isFlag(pointsAsNumbers),
            isFlag(label),
            isNumber(labelSize),
            isFlag(dark),
            isFlag(legend)
        )
        labels <- matchLabels(
            labels = labels,
            choices = eval(formals()[["labels"]])
        )
        dl(c(
            "reduction" = reduction,
            "dims" = deparse(dims)
        ))
        ## Note that we're not slotting interesting groups back into object
        ## here because we're allowing visualization of cluster identity, which
        ## isn't sample level.
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
        blacklist <- c("interestingGroups", "origIdent", "sampleId")
        supported <- setdiff(supported, blacklist)
        if (!isSubset(interestingGroups, supported)) {
            setdiff <- setdiff(interestingGroups, supported)
            abort(sprintf(
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
                toInlineString(setdiff, n = 5L),
                printString(supported)
            ))
        }
        data <- uniteInterestingGroups(
            object = data,
            interestingGroups = interestingGroups
        )
        ## Turn off labeling if there's only 1 cluster.
        if (hasLength(levels(data[["ident"]]), n = 1L)) {
            label <- FALSE
        }
        ## Set the x- and y-axis labels (e.g. t_SNE1, t_SNE2). We're setting
        ## this up internally as the first two columns in the data frame.
        axes <- colnames(data)[seq_len(2L)]
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym("x"),
                y = !!sym("y"),
                color = !!sym("interestingGroups")
            )
        )
        ## Points as numbers.
        if (isTRUE(pointsAsNumbers)) {
            ## Increase the size, if necessary.
            if (pointSize < 4L) {
                alertWarning("Increase pointSize to 4.")
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
        ## Label.
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
        ## Color.
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
        ## Improve the axis breaks.
        p <- p +
            scale_x_continuous(breaks = pretty_breaks(n = 4L)) +
            scale_y_continuous(breaks = pretty_breaks(n = 4L))
        ## Labels.
        if (is.list(labels)) {
            labels[["x"]] <- makeLabel(axes[[1L]])
            labels[["y"]] <- makeLabel(axes[[2L]])
            labels[["color"]] <- paste(interestingGroups, collapse = ":\n")
            labels[["fill"]] <- labels[["color"]]
            p <- p + do.call(what = labs, args = labels)
        }
        p
    }

args <- c(
    "dark",
    "dims",
    "label",
    "labelSize",
    "legend",
    "pointAlpha",
    "pointSize",
    "pointsAsNumbers",
    "reduction"
)
args1 <- c(args, "color")
args2 <- c(args, "discreteColor")
formals(`plotReducedDim,SingleCellExperiment`)[args1] <- .formalsList[args2]
rm(args, args1, args2)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("SingleCellExperiment"),
    definition = `plotReducedDim,SingleCellExperiment`
)



## Updated 2021-03-02.
`plotReducedDim,Seurat` <-  # nolint
    function(object, ...) {
        validObject(object)
        idents <- .seuratWhichIdents(object)
        dl(c("idents" = idents))
        plotReducedDim(object = as(object, "SingleCellExperiment"), ...)
    }



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("Seurat"),
    definition = `plotReducedDim,Seurat`
)



## Updated 2020-02-21.
`plotPCA,SingleCellExperiment` <-  # nolint
    function(object, ...) {
        plotReducedDim(object = object, reduction = "PCA", ...)
    }



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



## Updated 2020-02-21.
`plotTSNE,SingleCellExperiment` <-  # nolint
    function(object, ...) {
        plotReducedDim(object = object, reduction = "TSNE", ...)
    }



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("SingleCellExperiment"),
    definition = `plotTSNE,SingleCellExperiment`
)



## Updated 2020-02-21.
`plotTSNE,Seurat` <-  # nolint
    `plotTSNE,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("Seurat"),
    definition = `plotTSNE,Seurat`
)



## Updated 2020-02-21.
`plotUMAP,SingleCellExperiment` <-  # nolint
    function(object, ...) {
        plotReducedDim(object = object, reduction = "UMAP", ...)
    }



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature("SingleCellExperiment"),
    definition = `plotUMAP,SingleCellExperiment`
)



## Updated 2020-02-21.
`plotUMAP,Seurat` <-  # nolint
    `plotUMAP,SingleCellExperiment`



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature("Seurat"),
    definition = `plotUMAP,Seurat`
)

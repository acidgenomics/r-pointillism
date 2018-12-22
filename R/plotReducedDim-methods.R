#' @name plotReducedDim
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#' @inherit bioverbs::plotReducedDim
#' @inheritParams basejump::params
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
#'
#' We recommend installing this with [conda][]:
#'
#' ```
#' conda install -c conda-forge umap-learn
#' ```
#'
#' The UMAP module can be loaded in R using [reticulate][]. Reticulate works
#' reliably when setting `RETICULATE_PYTHON` to point to your conda python
#' binary. Export this variable in  `~/.Renviron`.
#'
#' See `base::Sys.getenv` for more information on the R system environment.
#'
#' [conda]: https://conda.io
#' [Python]: https://www.python.org
#' [UMAP]: https://github.com/lmcinnes/umap
#' [reticulate]: https://rstudio.github.io/reticulate/
#'
#' @seealso
#' - `Seurat::DimPlot`.
#' - [Seurat Mouse Cell Atlas vignette](https://satijalab.org/seurat/mca.html).
#'
#' @examples
#' data(seurat_small)
#' object <- seurat_small
#'
#' ## t-SNE
#' plotTSNE(object)
#' plotTSNE(object, pointsAsNumbers = TRUE, dark = TRUE, label = FALSE)
#'
#' ## UMAP
#' plotUMAP(object)
#'
#' ## PCA
#' plotPCA(object)
NULL



#' @importFrom bioverbs plotReducedDim
#' @aliases NULL
#' @export
bioverbs::plotReducedDim

#' @importFrom BiocGenerics plotPCA
#' @aliases NULL
#' @export
BiocGenerics::plotPCA

#' @importFrom bioverbs plotTSNE
#' @aliases NULL
#' @export
bioverbs::plotTSNE

#' @importFrom bioverbs plotUMAP
#' @aliases NULL
#' @export
bioverbs::plotUMAP



# Constructors =================================================================
plotReducedDim.SingleCellExperiment <-  # nolint
    function(
        object,
        reducedDim,
        dimsUse,
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
        assert(
            .hasIdent(object),
            isScalar(reducedDim),
            isIntegerish(dimsUse),
            hasLength(dimsUse, n = 2L),
            isString(interestingGroups, nullOK = TRUE),
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

        # Color by `ident` factor by default.
        if (is.null(interestingGroups)) {
            interestingGroups <- "ident"
        }

        data <- .fetchReducedDimData(
            object = object,
            reducedDim = reducedDim,
            dimsUse = dimsUse
        )
        assert(
            is(data, "DataFrame"),
            isSubset(
                x = c("x", "y", "centerX", "centerY", "interestingGroups"),
                y = colnames(data)
            )
        )

        # Color by `ident` factor by default (see above).
        if (interestingGroups == "ident") {
            data[["interestingGroups"]] <- data[["ident"]]
        }

        # Set the x- and y-axis labels (e.g. t_SNE1, t_SNE2).
        axes <- colnames(reducedDims(object)[[reducedDim]])[dimsUse]
        assert(isSubset(axes, colnames(data)))

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
            # Increase the size, if necessary.
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

        # Dark mode.
        if (isTRUE(dark)) {
            p <- p + theme_midnight()
        }

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        # Improve the axis breaks.
        p <- p +
            scale_x_continuous(breaks = pretty_breaks(n = 4L)) +
            scale_y_continuous(breaks = pretty_breaks(n = 4L))

        p
    }

formals(plotReducedDim.SingleCellExperiment)[c(
    "color",
    "dark",
    "dimsUse",
    "label",
    "labelSize",
    "legend",
    "pointAlpha",
    "pointSize",
    "pointsAsNumbers"
)] <- list(
    color = discreteColor,
    dark = dark,
    dimsUse = dimsUse,
    label = label,
    labelSize = labelSize,
    legend = legend,
    pointAlpha = pointAlpha,
    pointSize = pointSize,
    pointsAsNumbers = pointsAsNumbers
)



plotPCA.SingleCellExperiment <-  # nolint
    function() {
        do.call(
            what = plotReducedDim,
            args = matchArgsToDoCall(args = list(reducedDim = "PCA"))
        )

    }




plotTSNE.SingleCellExperiment <-  # nolint
    function() {
        do.call(
            what = plotReducedDim,
            args = matchArgsToDoCall(args = list(reducedDim = "TSNE"))
        )
    }



plotUMAP.SingleCellExperiment <-  # nolint
    function() {
        do.call(
            what = plotReducedDim,
            args = matchArgsToDoCall(args = list(reducedDim = "UMAP"))
        )
    }



# Formals ======================================================================
# Set the formals for the convenience functions.
f <- formals(plotReducedDim.SingleCellExperiment)
f <- f[setdiff(names(f), "reducedDim")]
formals(plotPCA.SingleCellExperiment) <- f
formals(plotTSNE.SingleCellExperiment) <- f
formals(plotUMAP.SingleCellExperiment) <- f
rm(f)



# Methods ======================================================================
#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("SingleCellExperiment"),
    definition = plotReducedDim.SingleCellExperiment
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("seurat"),
    definition = plotReducedDim.SingleCellExperiment
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("SingleCellExperiment"),
    definition = plotTSNE.SingleCellExperiment
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("seurat"),
    definition = plotTSNE.SingleCellExperiment
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature("SingleCellExperiment"),
    definition = plotUMAP.SingleCellExperiment
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature("seurat"),
    definition = plotUMAP.SingleCellExperiment
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("SingleCellExperiment"),
    definition = plotPCA.SingleCellExperiment
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("seurat"),
    definition = plotPCA.SingleCellExperiment
)

#' Plot Reduced Dimensions
#'
#' - t-SNE: **t**-distributed **S**tochastic **N**eighbor **E**mbedding.
#' - PCA: **P**rincipal **C**omponent **A**nalysis.
#' - UMAP: **U**niform **M**anifold **A**pproximation and **P**rojection.
#'
#' @name plotReducedDim
#' @family Plot Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom BiocGenerics plotPCA
#'
#' @inheritParams general
#'
#' @note [plotUMAP()] requires the Python dependency `umap-learn`. We recommend
#' installing this with conda: `conda install -c conda-forge umap-learn`.
#'
#' @seealso
#' - [Seurat::DimPlot()].
#' - [UMAP GitHub repo](https://github.com/lmcinnes/umap).
#' - [Seurat Mouse Cell Atlas vignette](https://satijalab.org/seurat/mca.html).
#'
#' @return `ggplot`.
#'
#' @examples
#' object <- sce_small
#'
#' # t-SNE
#' plotTSNE(object)
#' plotTSNE(object, pointsAsNumbers = TRUE, dark = TRUE, label = FALSE)
#'
#' # UMAP
#' plotUMAP(object)
#'
#' # PCA
#' plotPCA(object)
NULL



# Constructors =================================================================
.plotReducedDim <- function(
    object,
    reducedDim,
    dimsUse = c(1L, 2L),
    interestingGroups = "ident",
    color = getOption("pointillism.discrete.color", NULL),
    pointSize = getOption("pointillism.pointSize", 0.75),
    pointAlpha = getOption("pointillism.pointAlpha", 0.75),
    pointsAsNumbers = getOption("pointillism.pointsAsNumbers", FALSE),
    label = getOption("pointillism.label", TRUE),
    labelSize = getOption("pointillism.labelSize", 6L),
    dark = getOption("pointillism.dark", FALSE),
    legend = getOption("pointillism.legend", TRUE),
    title = NULL
) {
    .assertHasIdent(object)
    assert_is_a_string(reducedDim)
    assertIsImplicitInteger(dimsUse)
    assert_is_of_length(dimsUse, 2L)
    assert_is_a_string(interestingGroups)
    interestingGroups(object) <- interestingGroups
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_number(pointSize)
    assert_is_a_number(pointAlpha)
    assert_is_a_bool(pointsAsNumbers)
    assert_is_a_bool(label)
    assert_is_a_number(labelSize)
    assert_is_a_bool(dark)
    assert_is_a_bool(legend)
    assertIsAStringOrNULL(title)

    data <- .fetchReducedDimData(
        object = object,
        reducedDim = reducedDim,
        dimsUse = dimsUse
    )
    assert_is_all_of(data, "DataFrame")

    assert_is_subset(
        x = c("x", "y", "centerX", "centerY"),
        y = colnames(data)
    )

    # Set the x- and y-axis labels (e.g. tSNE1, tSNE2).
    axes <- camel(colnames(reducedDims(object)[[reducedDim]])[dimsUse])
    assert_is_subset(axes, colnames(data))

    p <- ggplot(
        data = as.data.frame(data),
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
        if (pointSize < 4L) pointSize <- 4L
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



.plotPCA <- function() {
    args <- as.list(match.call())[-1L]
    args[["reducedDim"]] <- "PCA"
    do.call(
        what = plotReducedDim,
        args = args
    )

}




.plotTSNE <- function() {
    args <- as.list(match.call())[-1L]
    args[["reducedDim"]] <- "UMAP"
    do.call(
        what = plotReducedDim,
        args = args
    )
}



.plotUMAP <- function() {
    args <- as.list(match.call())[-1L]
    args[["reducedDim"]] <- "UMAP"
    do.call(
        what = plotReducedDim,
        args = args
    )
}



# Formals ======================================================================
# Set the formals for the convenience functions.
f <- formals(.plotReducedDim)
f <- as.pairlist(f[setdiff(names(f), "reducedDim")])
formals(.plotPCA) <- f
formals(.plotTSNE) <- f
formals(.plotUMAP) <- f
rm(f)



# Methods ======================================================================
#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("SingleCellExperiment"),
    definition = .plotReducedDim
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature("seurat"),
    definition = getMethod("plotReducedDim", "SingleCellExperiment")
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("SingleCellExperiment"),
    definition = .plotTSNE
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature("seurat"),
    definition = getMethod("plotTSNE", "SingleCellExperiment")
)



#' @rdname plotReducedDim
#' @export
setMethod(
    "plotUMAP",
    signature("SingleCellExperiment"),
    .plotUMAP
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature("seurat"),
    definition = getMethod("plotUMAP", "SingleCellExperiment")
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("SingleCellExperiment"),
    definition = .plotPCA
)



#' @rdname plotReducedDim
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("seurat"),
    definition = getMethod("plotPCA", "SingleCellExperiment")
)

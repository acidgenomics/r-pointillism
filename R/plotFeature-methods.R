# TODO Add pointsAsNumbers support.
# TODO Is there a way to facet wrap these instead of using plot grid? Then
# we can easily support a title.
# FIXME We're using continuous color here, so the formal won't match...
# FIXME argument "color" is missing, with no default



#' @name plotFeature
#' @inherit bioverbs::plotFeature
#' @inheritParams basejump::params
#'
#' @param features `character`. Features to plot (e.g. gene expression, PC
#'   scores, number of genes detected).
#'
#' @seealso `Seurat::FeaturePlot()`.
#'
#' @return `ggplot` or `list`.
#'
#' @examples
#' data(seurat_small)
#' plotFeature(seurat_small, features = c("nUMI", "nGene"))
NULL



#' @importFrom bioverbs plotFeature
#' @aliases NULL
#' @export
bioverbs::plotFeature



# Note: We aren't using `title` argument here.
plotFeature.SingleCellExperiment <- function(
    object,
    features,
    reducedDim,
    color,
    pointSize,
    pointAlpha,
    pointsAsNumbers,
    label,
    labelSize,
    dark,
    legend
) {
    # Legacy arguments ---------------------------------------------------------
    # color
    if (identical(color, "auto")) {
        warning("Use `NULL` instead of `\"auto\"` for `color`")
        color <- NULL
    }

    # Assert checks ------------------------------------------------------------
    assert_is_character(features)
    # Sanitize input to camel case.
    features <- camel(features)
    reducedDim <- match.arg(reducedDim)
    assertIsColorScaleContinuousOrNULL(color)
    assert_is_a_number(pointSize)
    assert_is_a_number(pointAlpha)
    assert_is_a_bool(pointsAsNumbers)
    assert_is_a_bool(label)
    assert_is_a_number(labelSize)
    assert_is_a_bool(dark)
    assert_is_a_bool(legend)

    # Dark mode setup.
    if (isTRUE(dark)) {
        fill <- "black"
    } else {
        fill <- "white"
    }

    data <- .fetchReducedDimData(
        object = object,
        reducedDim = reducedDim
    )

    # Label the axes.
    axes <- colnames(data)[seq_len(2L)]

    # If the features are not defined, attempt to merge all reduced dims
    # information before stopping
    if (!all(features %in% colnames(data))) {
        reducedDimsData <- do.call(
            what = cbind,
            args = reducedDims(object)
        )
        assert_that(identical(rownames(data), rownames(reducedDimsData)))
        reducedDimsData <- camel(reducedDimsData)
        data <- data %>%
            .[, setdiff(colnames(.), colnames(reducedDimsData))] %>%
            cbind(reducedDimsData)
    }
    assert_is_subset(features, colnames(data))

    # FIXME Need to add pointsAsNumbers support.
    if (isTRUE(pointsAsNumbers)) {
        stop("pointsAsNumbers isn't supported yet")
    }

    plotlist <- lapply(features, function(feature) {
        p <- ggplot(
            data = as_tibble(data),
            mapping = aes(
                x = !!sym("x"),
                y = !!sym("y"),
                color = !!sym(feature)
            )
        ) +
            geom_point(
                alpha = pointAlpha,
                size = pointSize
            ) +
            labs(
                x = axes[[1L]],
                y = axes[[2L]],
                title = feature
            )

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

        if (isTRUE(dark)) {
            p <- p + theme_midnight()
            if (is.null(color)) {
                color <- darkMarkerColors
            }
        }

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        if (!isTRUE(legend)) {
            p <- p + guides(color = "none")
        }

        p
    })

    # Return ---------------------------------------------------------------
    if (length(features) > 1L) {
        plot_grid(plotlist = plotlist) +
            theme(
                plot.background = element_rect(color = NA, fill = fill)
            )
    } else {
        plotlist[[1L]]
    }
}

formals(plotFeature.SingleCellExperiment)[c(
    "color",
    "dark",
    "expression",
    "label",
    "labelSize",
    "legend",
    "pointAlpha",
    "pointSize",
    "pointsAsNumbers",
    "reducedDim"
)] <- list(
    color = continuousColor,
    dark = dark,
    expression = expression,
    label = label,
    labelSize = labelSize,
    legend = legend,
    pointAlpha = pointAlpha,
    pointSize = pointSize,
    pointsAsNumbers = pointsAsNumbers,
    reducedDim = reducedDim
)




#' @rdname plotFeature
#' @export
setMethod(
    f = "plotFeature",
    signature = signature("SingleCellExperiment"),
    definition = plotFeature.SingleCellExperiment
)



#' @rdname plotFeature
#' @export
setMethod(
    f = "plotFeature",
    signature = signature("seurat"),
    definition = plotFeature.SingleCellExperiment
)

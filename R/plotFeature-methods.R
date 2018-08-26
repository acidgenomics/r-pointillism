#' Plot Feature
#'
#' @name plotFeature
#' @family Plot Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param features `character`. Features to plot (e.g. gene expression, PC
#'   scores, number of genes detected).
#'
#' @seealso [Seurat::FeaturePlot()].
#'
#' @return `ggplot` or `list`.
#'
#' @examples
#' plotFeature(sce_small, features = c("nUMI", "nGene", "mitoRatio"))
NULL



#' @rdname plotFeature
#' @export
setMethod(
    "plotFeature",
    signature("SingleCellExperiment"),
    function(
        object,
        features,
        reducedDim = c("TSNE", "UMAP"),
        color = getOption("pointillism.discrete.color", NULL),
        pointSize = getOption("pointillism.pointSize", 0.75),
        pointAlpha = getOption("pointillism.pointAlpha", 0.75),
        label = getOption("pointillism.label", TRUE),
        labelSize = getOption("pointillism.labelSize", 6L),
        dark = getOption("pointillism.dark", FALSE),
        legend = getOption("pointillism.legend", TRUE)
    ) {
        assert_is_character(features)
        # Sanitize input to camel case.
        features <- camel(features)
        # Legacy support for `color = "auto"`
        if (identical(color, "auto")) {
            warning("Use `NULL` instead of `\"auto\"` for `color`")
            color <- NULL
        }
        assertIsColorScaleContinuousOrNULL(color)
        reducedDim <- match.arg(reducedDim)

        if (isTRUE(dark)) {
            fill <- "black"
        } else {
            fill <- "white"
        }

        data <- .fetchReducedDimData(object, reducedDim = reducedDim)
        axes <- colnames(data)[seq_len(2L)]

        # If the features are not defined, attempt to merge all reduced dims
        # information before stopping
        if (!all(features %in% colnames(data))) {
            reducedDimsData <- do.call(
                what = cbind,
                args = reducedDims(object)
            )
            stopifnot(identical(rownames(data), rownames(reducedDimsData)))
            reducedDimsData <- camel(reducedDimsData)
            data <- data %>%
                .[, setdiff(colnames(.), colnames(reducedDimsData))] %>%
                cbind(reducedDimsData)
        }
        assert_is_subset(features, colnames(data))

        plotlist <- lapply(features, function(feature) {
            p <- ggplot(
                data = data,
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
                    plot.background = element_rect(
                        color = NA,
                        fill = fill
                    )
                )
        } else {
            plotlist[[1L]]
        }
    }
)



#' @rdname plotFeature
#' @export
setMethod(
    "plotFeature",
    signature("seurat"),
    getMethod("plotFeature", "SingleCellExperiment")
)

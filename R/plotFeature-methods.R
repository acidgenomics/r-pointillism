## FIXME Fix handling for reducedDim matrices without names (e.g. UMAP) for
## monocle3



#' @name plotFeature
#' @inherit bioverbs::plotFeature
#' @note Updated 2019-07-31.
#'
#' @inheritParams acidplots::params
#' @inheritParams acidroxygen::params
#' @param features `character`. Features to plot (e.g. gene expression, PC
#'   scores, number of genes detected).
#' @param ... Additional arguments.
#'
#' @seealso [Seurat::FeaturePlot()].
#'
#' @return `ggplot` (1 feature) or `list` (multiple features).
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' plotFeature(object, features = c("nCount_RNA", "nFeature_RNA"))
NULL



#' @rdname plotFeature
#' @name plotFeature
#' @importFrom bioverbs plotFeature
#' @usage plotFeature(object, ...)
#' @export
NULL



## Note: We aren't using `title` argument here.
## Updated 2019-07-31.
`plotFeature,SingleCellExperiment` <-  # nolint
    function(
        object,
        features,
        reduction,
        color,
        pointSize,
        pointAlpha,
        pointsAsNumbers,
        label,
        labelSize,
        dark,
        legend
    ) {
        ## Legacy arguments ----------------------------------------------------
        ## color
        if (identical(color, "auto")) {
            warning("Use `NULL` instead of `\"auto\"` for `color`")
            color <- NULL
        }

        ## Assert checks -------------------------------------------------------
        assert(isCharacter(features))
        reduction <- match.arg(reduction)
        assert(
            isGGScale(
                x = color,
                scale = "continuous",
                aes = "colour",
                nullOK = TRUE
            ),
            isNumber(pointSize),
            isNumber(pointAlpha),
            isFlag(pointsAsNumbers),
            isFlag(label),
            isNumber(labelSize),
            isFlag(dark),
            isFlag(legend)
        )

        ## Dark mode setup.
        if (isTRUE(dark)) {
            fill <- "black"
        } else {
            fill <- "white"
        }

        data <- .fetchReductionData(
            object = object,
            reduction = reduction
        )

        ## Label the axes.
        axes <- colnames(data)[seq_len(2L)]

        ## If the features are not defined, attempt to merge all reduced dims
        ## information before stopping.
        if (!all(features %in% colnames(data))) {
            reductionData <- do.call(
                what = cbind,
                args = reducedDims(object)
            )
            assert(identical(rownames(data), rownames(reductionData)))
            ## FIXME We need to rethink this approach for monocle3.
            reductionData <- camel(reductionData)
            data <- data %>%
                .[, setdiff(colnames(.), colnames(reductionData))] %>%
                cbind(reductionData)
        }
        assert(isSubset(features, colnames(data)))

        ## Need to add pointsAsNumbers support.
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
                p <- p + acid_theme_dark()
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

        ## Return --------------------------------------------------------------
        if (length(features) > 1L) {
            plot_grid(plotlist = plotlist) +
                theme(
                    plot.background = element_rect(color = NA, fill = fill)
                )
        } else {
            plotlist[[1L]]
        }
    }

formals(`plotFeature,SingleCellExperiment`)[c(
    "color",
    "dark",
    "expression",
    "label",
    "labelSize",
    "legend",
    "pointAlpha",
    "pointSize",
    "pointsAsNumbers",
    "reduction"
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
    reduction = reduction
)



#' @rdname plotFeature
#' @export
setMethod(
    f = "plotFeature",
    signature = signature("SingleCellExperiment"),
    definition = `plotFeature,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotFeature,Seurat` <-  # nolint
    `plotFeature,SingleCellExperiment`



#' @rdname plotFeature
#' @export
setMethod(
    f = "plotFeature",
    signature = signature("Seurat"),
    definition = `plotFeature,Seurat`
)

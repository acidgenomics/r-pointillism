# TODO Consolidate these defaults with `.plotReducedDim`.
# FIXME Keep the original gene input, to use for the labels.
# Need to remap an additional remapping function here I think.



#' Plot Cell-Type-Specific Gene Markers
#'
#' Visualize gene markers on a reduced dimension plot (e.g. t-SNE, UMAP).
#'
#' @name plotMarker
#' @family Plot Functions
#' @author Michael Steinbaugh, Rory Kirchner
#' @export
#'
#' @inheritParams general
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' title <- "most abundant genes"
#' genes <- counts(sce_small) %>%
#'     Matrix::rowSums(.) %>%
#'     sort(decreasing = TRUE) %>%
#'     head(n = 4L) %>%
#'     names()
#' glimpse(genes)
#' plotMarker(sce_small, genes = genes[[1L]])
#' plotMarker(
#'     object = sce_small,
#'     genes = genes,
#'     expression = "mean",
#'     pointsAsNumbers = TRUE,
#'     dark = TRUE,
#'     label = FALSE,
#'     title = title
#' )
NULL



.plotMarker <- function(
    object,
    genes,
    reducedDim = c("TSNE", "UMAP"),
    expression = c("mean", "sum"),
    color = getOption(
        "pointillism.continuous.color",
        ggplot2::scale_color_gradient(
            low = "orange",
            high = "purple"
        )
    ),
    pointSize = getOption("pointillism.pointSize", 0.75),
    pointAlpha = getOption("pointillism.pointAlpha", 0.8),
    pointsAsNumbers = getOption("pointillism.pointsAsNumbers", FALSE),
    label = getOption("pointillism.label", TRUE),
    labelSize = getOption("pointillism.labelSize", 6L),
    dark = getOption("pointillism.dark", FALSE),
    legend = getOption("pointillism.legend", TRUE),
    title = TRUE
) {
    # Legacy arguments -----------------------------------------------------
    # color
    if (identical(color, "auto")) {
        message("Use `color = NULL` instead of `auto`")
        color <- NULL
    }

    # Assert checks --------------------------------------------------------
    validObject(object)
    assert_is_character(genes)
    geneNames <- mapGenesToSymbols(object, genes)
    reducedDim <- match.arg(reducedDim)
    expression <- match.arg(expression)
    assertIsColorScaleContinuousOrNULL(color)
    assert_is_a_number(pointSize)
    assert_is_a_number(pointAlpha)
    assert_is_a_bool(pointsAsNumbers)
    assert_is_a_bool(label)
    assert_is_a_number(labelSize)
    assert_is_a_bool(dark)
    assert_is_a_bool(legend)
    assert_is_any_of(title, c("character", "logical", "NULL"))
    if (is.character(title)) {
        assert_is_a_string(title)
    }

    # Fetch reduced dimension data
    data <- .fetchReducedDimExpressionData(
        object = object,
        genes = genes,
        reducedDim = reducedDim
    )
    assert_is_all_of(data, "DataFrame")

    # Get the axis labels.
    axes <- colnames(data)[seq_len(2L)]
    stopifnot(all(grepl("\\d+$", axes)))

    requiredCols <- c(
        axes,
        "centerX",
        "centerY",
        "ident",
        "mean",
        "sum",
        "x",
        "y"
    )
    assert_is_subset(requiredCols, colnames(data))

    p <- ggplot(
        data = as.data.frame(data),
        mapping = aes(
            x = !!sym("x"),
            y = !!sym("y"),
            color = !!sym(expression)
        )
    )

    # Titles
    subtitle <- NULL
    if (isTRUE(title)) {
        if (is_a_string(geneNames)) {
            title <- geneNames
        } else {
            title <- NULL
            subtitle <- geneNames
            # Limit to the first 5 markers
            if (length(subtitle) > 5L) {
                subtitle <- c(subtitle[1L:5L], "...")
            }
            subtitle <- toString(subtitle)
        }
    } else if (identical(title, FALSE)) {
        title <- NULL
    }
    p <- p +
        labs(
            x = axes[[1L]],
            y = axes[[2L]],
            title = title,
            subtitle = subtitle
        )

    # Customize legend.
    if (isTRUE(legend)) {
        if (is_a_string(genes)) {
            guideTitle <- "logcounts"
        } else {
            guideTitle <- paste0(
                "logcounts", "\n",
                "(", expression, ")"
            )
        }
        p <- p + guides(color = guide_colorbar(title = guideTitle))
    } else {
        p <- p + guides(color = "none")
    }

    if (isTRUE(pointsAsNumbers)) {
        if (pointSize < 4L) pointSize <- 4L
        p <- p +
            geom_text(
                mapping = aes(
                    x = !!sym("x"),
                    y = !!sym("y"),
                    label = !!sym("ident"),
                    color = !!sym(expression)
                ),
                alpha = pointAlpha,
                size = pointSize
            )
    } else {
        p <- p +
            geom_point(
                alpha = pointAlpha,
                size = pointSize
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
        if (is.null(color)) {
            color <- darkMarkerColors
        }
    }

    if (is(color, "ScaleContinuous")) {
        p <- p + color
    }

    p
}



#' @rdname plotMarker
#' @export
setMethod(
    "plotMarker",
    signature("SingleCellExperiment"),
    .plotMarker
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotMarker",
    signature("seurat"),
    getMethod("plotMarker", "SingleCellExperiment")
)



# Features =====================================================================
# TODO Add pointsAsNumbers support.
# TODO Is there a way to facet wrap these instead of using plot grid? Then
# we can easily support a title.


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



# FIXME Simpilfy the formals here to match `plotReducedDim`.
# NOTE We aren't using `title` argument here.
# FIXME We're using continuous color here, so the formal won't match...
.plotFeature <- function(
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
        stopifnot(identical(rownames(data), rownames(reducedDimsData)))
        reducedDimsData <- camel(reducedDimsData)
        data <- data %>%
            .[, setdiff(colnames(.), colnames(reducedDimsData))] %>%
            cbind(reducedDimsData)
    }
    assert_is_subset(features, colnames(data))

    # FIXME Need to add pointsAsNumbers support here
    if (isTRUE(pointsAsNumbers)) {
        stop("pointsAsNumbers isn't supported yet")
    }

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
                plot.background = element_rect(color = NA, fill = fill)
            )
    } else {
        plotlist[[1L]]
    }
}




#' @rdname plotFeature
#' @export
setMethod(
    f = "plotFeature",
    signature = signature("SingleCellExperiment"),
    definition = .plotFeature
)



#' @rdname plotFeature
#' @export
setMethod(
    f = "plotFeature",
    signature = signature("seurat"),
    definition = getMethod(
        f = "plotFeature",
        signature = signature("SingleCellExperiment")
    )
)

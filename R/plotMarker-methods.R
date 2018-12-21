#' @name plotMarker
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit bioverbs::plotMarker
#' @inheritParams basejump::params
#' @examples
#' data(seurat_small)
#' object <- seurat_small
#' title <- "most abundant genes"
#' genes <- counts(object) %>%
#'     Matrix::rowSums(.) %>%
#'     sort(decreasing = TRUE) %>%
#'     head(n = 4L) %>%
#'     names()
#' str(genes)
#'
#' ## Default appearance.
#' plotMarker(object, genes = genes[[1L]])
#'
#' ## Dark mode with viridis palette.
#' plotMarker(
#'     object = object,
#'     genes = genes,
#'     expression = "mean",
#'     pointsAsNumbers = TRUE,
#'     color = ggplot2::scale_color_viridis_c(),
#'     dark = TRUE,
#'     label = FALSE,
#'     title = title
#' )
NULL



#' @importFrom bioverbs plotMarker
#' @aliases NULL
#' @export
bioverbs::plotMarker



plotMarker.SingleCellExperiment <- function(
    object,
    genes,
    reducedDim,
    expression,
    color,
    pointSize,
    pointAlpha,
    pointsAsNumbers,
    label,
    labelSize,
    dark,
    legend,
    title = TRUE
) {
    # Legacy arguments -----------------------------------------------------
    # color
    if (identical(color, "auto")) {
        stop("Use `color = NULL` instead of `auto`.")
    }

    # Assert checks --------------------------------------------------------
    object <- as(object, "SingleCellExperiment")
    assert_is_character(genes)
    # FIXME This step is breaking.
    geneNames <- mapGenesToSymbols(object, genes)
    assert_is_scalar(reducedDim)
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
    assert_that(all(grepl("\\d+$", axes)))

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
        data = as_tibble(data),
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
formals(plotMarker.SingleCellExperiment)[c(
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



#' @rdname plotMarker
#' @export
setMethod(
    f = "plotMarker",
    signature = signature("SingleCellExperiment"),
    definition = plotMarker.SingleCellExperiment
)



#' @rdname plotMarker
#' @export
setMethod(
    f = "plotMarker",
    signature = signature("seurat"),
    definition = plotMarker.SingleCellExperiment
)

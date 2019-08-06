#' @name plotMarker
#' @author Michael Steinbaugh, Rory Kirchner
#' @note Updated 2019-08-02.
#'
#' @inherit bioverbs::plotMarker
#'
#' @inheritParams acidplots::params
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
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
#' genes <- counts(object) %>%
#'     Matrix::rowSums(.) %>%
#'     sort(decreasing = TRUE) %>%
#'     head(n = 4L) %>%
#'     names()
#' print(genes)
#' plotMarker(object, genes = genes)
#'
#' ## cell_data_set ====
#' object = cell_data_set
#' genes <- counts(object) %>%
#'     Matrix::rowSums(.) %>%
#'     sort(decreasing = TRUE) %>%
#'     head(n = 4L) %>%
#'     names()
#' print(genes)
#' plotMarker(object, genes = genes)
NULL



#' @rdname plotMarker
#' @name plotMarker
#' @importFrom bioverbs plotMarker
#' @usage plotMarker(object, ...)
#' @export
NULL



## Updated 2019-08-02.
`plotMarker,SingleCellExperiment` <-  # nolint
    function(
        object,
        genes,
        reduction,
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
        ## Legacy arguments ----------------------------------------------------
        ## color
        if (identical(color, "auto")) {
            stop("Use `color = NULL` instead of `auto`.")
        }

        ## Assert checks -------------------------------------------------------
        assert(
            isCharacter(genes),
            isScalar(reduction),
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
            isFlag(legend),
            isAny(title, c("character", "logical", "NULL"))
        )
        geneNames <- mapGenesToSymbols(object, genes)
        expression <- match.arg(expression)
        if (is.character(title)) {
            assert(isString(title))
        }

        ## Fetch reduced dimension data.
        data <- .fetchReductionExpressionData(
            object = object,
            genes = genes,
            reduction = reduction
        )
        assert(is(data, "DataFrame"))

        ## Get the axis labels.
        axes <- colnames(data)[seq_len(2L)]
        assert(all(grepl("\\d+$", axes)))

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
        assert(isSubset(requiredCols, colnames(data)))

        p <- ggplot(
            data = as_tibble(data),
            mapping = aes(
                x = !!sym("x"),
                y = !!sym("y"),
                color = !!sym(expression)
            )
        )

        ## Titles
        subtitle <- NULL
        if (isTRUE(title)) {
            if (isString(geneNames)) {
                title <- geneNames
            } else {
                title <- NULL
                subtitle <- geneNames
                ## Limit to the first 5 markers
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

        ## Customize legend.
        if (isTRUE(legend)) {
            if (isString(genes)) {
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

        ## Dark mode.
        if (isTRUE(dark)) {
            p <- p + acid_theme_dark()
            if (is.null(color)) {
                color <- darkMarkerColors
            }
        }

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        p
    }

formals(`plotMarker,SingleCellExperiment`)[c(
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



#' @rdname plotMarker
#' @export
setMethod(
    f = "plotMarker",
    signature = signature("SingleCellExperiment"),
    definition = `plotMarker,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotMarker,Seurat` <-  # nolint
    `plotMarker,SingleCellExperiment`



#' @rdname plotMarker
#' @export
setMethod(
    f = "plotMarker",
    signature = signature("Seurat"),
    definition = `plotMarker,Seurat`
)

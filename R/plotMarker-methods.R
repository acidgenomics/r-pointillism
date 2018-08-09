#' Plot Cell-Type-Specific Gene Markers
#'
#' Visualize gene markers on a reduced dimension plot (e.g. t-SNE, UMAP).
#'
#' @name plotMarker
#' @family Clustering Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' # MHC class II genes
#' title <- "MHC class II"
#' genes <- rownames(sce_small)[which(grepl(
#'     pattern = "^major histocompatibility complex, class II",
#'     x = rowData(sce_small)$description
#' ))]
#' print(genes)
#'
#' plotMarker(sce_small, genes = genes)
#' plotMarker(
#'     object = sce_small,
#'     genes = genes,
#'     expression = "sum",
#'     pointsAsNumbers = TRUE,
#'     dark = TRUE,
#'     label = FALSE,
#'     title = title
#' )
NULL



# Constructors =================================================================
# Strip everything except the x-axis text labels
.minimalAxis <- function() {
    theme(
        axis.line = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        title = element_blank()
    )
}



# Methods ======================================================================
#' @rdname plotMarker
#' @export
setMethod(
    "plotMarker",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        reducedDim = c("TSNE", "UMAP", "PCA"),
        expression = c("mean", "median", "sum"),
        color = getOption("bcbio.discrete.color", NULL),
        pointSize = getOption("bcbio.pointSize", 0.75),
        pointAlpha = getOption("bcbio.pointAlpha", 0.8),
        pointsAsNumbers = FALSE,
        label = getOption("bcbio.label", TRUE),
        labelSize = getOption("bcbio.labelSize", 6L),
        dark = getOption("bcbio.dark", FALSE),
        grid = getOption("bcbio.grid", FALSE),
        legend = getOption("bcbio.legend", TRUE),
        aspectRatio = getOption("bcbio.aspectRatio", 1L),
        title = TRUE
    ) {
        assert_is_character(genes)
        assert_has_no_duplicates(genes)
        assert_is_subset(genes, rownames(object))
        reducedDim <- match.arg(reducedDim)
        expression <- match.arg(expression)
        # Legacy support for `color = "auto"`
        if (identical(color, "auto")) {
            color <- NULL
        }
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
        axes <- colnames(data)[seq_len(2L)]

        if (isTRUE(.useGene2symbol(object))) {
            g2s <- gene2symbol(object)
            symbols <- g2s[
                match(genes, g2s[["geneID"]]),
                "geneName",
                drop = TRUE
            ]
            stopifnot(!any(is.na(symbols)))
            genes <- make.unique(symbols)
        }
        genes <- sort(unique(genes))

        requiredCols <- c(
            axes,
            "x",
            "y",
            "centerX",
            "centerY",
            "mean",
            "median",
            "ident",
            "sum"
        )
        assert_is_subset(requiredCols, colnames(data))

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("x"),
                y = !!sym("y"),
                color = !!sym(expression)
            )
        )

        # Titles
        subtitle <- NULL
        if (isTRUE(title)) {
            if (is_a_string(genes)) {
                title <- genes
            } else {
                title <- NULL
                subtitle <- genes
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

        # Customize legend
        if (isTRUE(legend)) {
            if (is_a_string(genes)) {
                guideTitle <- "expression"
            } else {
                guideTitle <- expression
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

        # Color palette
        if (isTRUE(dark)) {
            theme <- theme_midnight
            if (is.null(color)) {
                color <- darkMarkerColors
            }
        } else {
            theme <- theme_paperwhite
            if (is.null(color)) {
                color <- lightMarkerColors
            }
        }
        p <- p +
            theme(
                aspect_ratio = aspectRatio,
                grid = grid
            )

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        p
    }
)

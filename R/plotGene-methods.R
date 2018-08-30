#' Plot Gene
#'
#' Visualize genes on a dot or violin plot.
#'
#' @name plotGene
#' @family Plot Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump plotGene
#' @export
#'
#' @inheritParams general
#' @inheritParams ggplot2::geom_violin
#' @param colMin `scalar numeric`. Minimum scaled average expression threshold.
#'   Everything smaller will be set to this.
#' @param colMax `scalar numeric`. Maximum scaled average expression threshold.
#'   Everything larger will be set to this.
#' @param dotMin `scalar numeric`. The fraction of cells at which to draw the
#'   smallest dot. All cell groups with less than this expressing the given gene
#'   will have no dot drawn.
#' @param dotScale `scalar numeric`. Scale the size of the points, similar to
#'   `cex`.
#' @param geom `string`. Plot type. Uses [match.arg()] to pick the type.
#'   Currently supports "`dot`" and "`violin`".
#'
#' @seealso
#' - [Seurat::DotPlot()].
#' - [Seurat::VlnPlot()].
#' - [Seurat::RidgePlot()].
#'
#' @return `ggplot`.
#'
#' @examples
#' object <- sce_small
#'
#' # Plotting with either gene IDs or gene names (symbols) works.
#' geneIDs <- head(rownames(object), n = 4L)
#' print(geneIDs)
#' geneNames <- head(as.character(rowRanges(object)$geneName), n = 4L)
#' print(geneNames)
#'
#' # Per sample mode enabled.
#' plotDot(object, genes = geneNames)
#' plotViolin(object, genes = geneNames)
#'
#' # Per sample mode disabled.
#' plotDot(object, genes = geneIDs, perSample = FALSE)
#' plotViolin(object, genes = geneIDs, perSample = FALSE)
NULL



# plotGene =====================================================================
.plotGene <- function(
    object,
    genes,
    geom = c("dot", "violin"),
    perSample = TRUE,
    color,
    legend = getOption("pointillism.legend", TRUE),
    title = NULL
) {
    validObject(object)
    geom <- match.arg(geom)
    if (geom == "dot") {
        what <- plotDot
    } else if (geom == "violin") {
        what <- plotViolin
    }
    args <- as.list(sys.call(which = -1L))[-1L]
    args[["geom"]] <- NULL
    do.call(what = what, args = args)
}



# plotDot ======================================================================
#' Min Max
#' @seealso [Seurat:::MinMax()].
#' @noRd
.minMax <- function(data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    data2
}



#' Percent Above
#' @seealso [Seurat:::PercentAbove()].
#' @noRd
.percentAbove <- function(x, threshold) {
    length(x[x > threshold]) / length(x)
}



# FIXME Need to fix the color scale argument here. It's continuous.
# FIXME Share this default with `plotMarker`.
.plotDot <- function(
    object,
    genes,
    perSample,
    colMin = -2.5,
    colMax = 2.5,
    dotMin = 0L,
    dotScale = 6L,
    color = getOption(
        "pointillism.continuous.color",
        ggplot2::scale_color_gradient(
            low = "orange",
            high = "purple"
        )
    ),
    legend,
    title
) {
    validObject(object)
    .assertHasIdent(object)
    assert_is_character(genes)
    assert_is_a_bool(perSample)
    assert_is_a_number(colMin)
    assert_is_a_number(colMax)
    assert_is_a_number(dotMin)
    assert_is_a_number(dotScale)
    assertIsColorScaleContinuousOrNULL(color)
    assert_is_a_bool(legend)
    assertIsAStringOrNULL(title)

    # Fetch the gene expression data.
    data <- .fetchGeneData(
        object = object,
        genes = genes,
        assay = "logcounts",
        metadata = TRUE
    )

    # Prepare data for ggplot.
    cols <- c("geneName", "sampleName", "ident")
    data <- data %>%
        as("tbl_df") %>%
        group_by(!!!syms(cols)) %>%
        summarize(
            avgExp = mean(expm1(!!sym("logcounts"))),
            # Consider making threshold user definable.
            pctExp = .percentAbove(!!sym("logcounts"), threshold = 0L)
        ) %>%
        ungroup() %>%
        mutate(geneName = as.factor(!!sym("geneName"))) %>%
        group_by(!!sym("geneName")) %>%
        mutate(
            avgExpScale = scale(!!sym("avgExp")),
            avgExpScale = .minMax(
                !!sym("avgExpScale"),
                max = colMax,
                min = colMin
            )
        ) %>%
        arrange(!!!syms(cols), .by_group = TRUE)

    # Apply our `dotMin` threshold.
    data[["pctExp"]][data[["pctExp"]] < dotMin] <- NA

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym("geneName"),
            y = !!sym("ident")
        )
    ) +
        geom_point(
            mapping = aes(
                color = !!sym("avgExpScale"),
                size = !!sym("pctExp")
            ),
            show.legend = legend
        ) +
        scale_radius(range = c(0L, dotScale)) +
        labs(
            x = NULL,
            y = NULL
        )

    # Handling step for multiple samples, if desired.
    if (
        isTRUE(perSample) &&
        isTRUE(.hasMultipleSamples(object))
    ) {
        p <- p +
            facet_wrap(
                facets = vars(!!sym("sampleName"))
            )
    }

    if (is(color, "ScaleContinuous")) {
        p <- p + color
    }

    p
}



.plotViolin <- function(
    object,
    genes,
    perSample,
    scale = c("count", "width", "area"),
    color = getOption("pointillism.discrete.color", NULL),
    legend,
    title
) {
    validObject(object)
    assert_is_character(genes)
    assert_is_a_bool(perSample)
    scale <- match.arg(scale)
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_bool(legend)
    assertIsAStringOrNULL(title)

    # Fetch the gene expression data.
    data <- .fetchGeneData(
        object = object,
        genes = genes,
        assay = "logcounts",
        metadata = TRUE
    )

    # Handling step for multiple samples, if desired.
    if (
        isTRUE(perSample) &&
        isTRUE(.hasMultipleSamples(object))
    ) {
        x <- "sampleName"
        interestingGroups <- interestingGroups(object)
        if (
            is.null(interestingGroups) ||
            interestingGroups == "ident"
        ) {
            interestingGroups <- "sampleName"
        }
        colorMapping <- "interestingGroups"
        colorLabs <- paste(interestingGroups, collapse = ":\n")
    } else {
        x <- "ident"
        colorMapping <- x
        colorLabs <- x
    }

    p <- ggplot(
        data = as.data.frame(data),
        mapping = aes(
            x = !!sym(x),
            y = !!sym("logcounts"),
            color = !!sym(colorMapping)
        )
    ) +
        geom_jitter(
            show.legend = legend
        ) +
        geom_violin(
            fill = NA,
            scale = scale,
            adjust = 1L,
            show.legend = legend,
            trim = TRUE
        ) +
        # Note that `scales = free_y` will hide the x-axis for some plots.
        labs(
            title = title,
            color = colorLabs
        )

    # Handling step for multiple samples, if desired.
    if (
        isTRUE(perSample) &&
        isTRUE(.hasMultipleSamples(object))
    ) {
        p <- p +
            facet_grid(
                rows = vars(!!sym("ident")),
                cols = vars(!!sym("geneName")),
                scales = "free_y"
            )
    } else {
        p <- p +
            facet_wrap(
                facets = vars(!!sym("geneName")),
                scales = "free_y"
            )
    }

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    p
}



# Methods ======================================================================
#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("SingleCellExperiment"),
    .plotGene
)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("seurat"),
    getMethod("plotGene", "SingleCellExperiment")
)



#' @rdname plotGene
#' @export
setMethod(
    "plotDot",
    signature("SingleCellExperiment"),
    .plotDot
)



#' @rdname plotGene
#' @export
setMethod(
    "plotDot",
    signature("seurat"),
    getMethod("plotDot", "SingleCellExperiment")
)



#' @rdname plotGene
#' @export
setMethod(
    "plotViolin",
    signature("SingleCellExperiment"),
    .plotViolin
)



#' @rdname plotGene
#' @export
setMethod(
    "plotViolin",
    signature("seurat"),
    getMethod("plotViolin", "SingleCellExperiment")
)

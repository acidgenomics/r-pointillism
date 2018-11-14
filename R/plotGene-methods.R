#' Plot Gene
#'
#' Visualize genes on a dot or violin plot.
#'
#' @name plotGene
#'
#' @inheritParams basejump::params
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
#' @param geom `string`. Plot type. Uses [base::match.arg()] to pick the type.
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
#' data(seurat_small)
#' object <- seurat_small
#'
#' ## Plotting with either gene IDs or gene names (symbols) works.
#' geneIDs <- head(rownames(object), n = 4L)
#' print(geneIDs)
#' geneNames <- head(as.character(rowRanges(object)$geneName), n = 4L)
#' print(geneNames)
#'
#' ## Per sample mode enabled.
#' plotDot(object, genes = geneNames, perSample = TRUE)
#' plotViolin(object, genes = geneNames, perSample = TRUE)
#'
#' ## Per sample mode disabled.
#' plotDot(object, genes = geneIDs, perSample = FALSE)
#' plotViolin(object, genes = geneIDs, perSample = FALSE)
NULL



#' @importFrom basejump plotGene
#' @aliases NULL
#' @export
basejump::plotGene



# plotGene =====================================================================
plotGene.SingleCellExperiment <- function(
    object,
    genes,
    geom = c("dot", "violin"),
    perSample = TRUE,
    legend,
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
formals(plotGene.SingleCellExperiment)[["legend"]] <- legend



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



plotDot.SingleCellExperiment <- function(
    object,
    genes,
    perSample = TRUE,
    colMin = -2.5,
    colMax = 2.5,
    dotMin = 0L,
    dotScale = 6L,
    color,
    legend,
    title = NULL
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
    assertIsStringOrNULL(title)

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
        as_tibble() %>%
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
formals(plotDot.SingleCellExperiment)[c(
    "color",
    "legend"
)] <- list(
    color = continuousColorPurpleOrange,
    legend = legend
)



plotViolin.SingleCellExperiment <- function(
    object,
    genes,
    perSample = TRUE,
    scale = c("count", "width", "area"),
    color,
    legend,
    title = NULL
) {
    validObject(object)
    assert_is_character(genes)
    assert_is_a_bool(perSample)
    scale <- match.arg(scale)
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_bool(legend)
    assertIsStringOrNULL(title)

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
        data = as_tibble(data),
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
formals(plotViolin.SingleCellExperiment)[c(
    "color",
    "legend"
)] <- list(
    color = discreteColor,
    legend = legend
)



# Methods ======================================================================
#' @rdname plotGene
#' @export
setMethod(
    f = "plotGene",
    signature = signature("SingleCellExperiment"),
    definition = plotGene.SingleCellExperiment
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotGene",
    signature = signature("seurat"),
    getMethod(
        f ="plotGene",
        signature = signature("SingleCellExperiment")
    )
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotDot",
    signature = signature("SingleCellExperiment"),
    definition = plotDot.SingleCellExperiment
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotDot",
    signature = signature("seurat"),
    getMethod(
        f = "plotDot",
        signature("SingleCellExperiment")
    )
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotViolin",
    signature = signature("SingleCellExperiment"),
    definition = plotViolin.SingleCellExperiment
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotViolin",
    signature = signature("seurat"),
    getMethod(
        f = "plotViolin",
        signature("SingleCellExperiment")
    )
)

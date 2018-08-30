# FIXME Consolidate this into plotGene.



#' Plot Dot
#'
#' @name plotDot
#' @family Plot Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param colMin `scalar numeric`. Minimum scaled average expression threshold.
#'   Everything smaller will be set to this.
#' @param colMax `scalar numeric`. Maximum scaled average expression threshold.
#'   Everything larger will be set to this.
#' @param dotMin `scalar numeric`. The fraction of cells at which to draw the
#'   smallest dot. All cell groups with less than this expressing the given gene
#'   will have no dot drawn.
#' @param dotScale `scalar numeric`. Scale the size of the points, similar to
#'   `cex`.
#'
#' @seealso Modified version of [Seurat::DotPlot()].
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
#' plotDot(object, genes = geneIDs)
#' plotDot(object, genes = geneNames)
#'
#' # Per sample mode can be disabled.
#' plotDot(object, genes = geneIDs, perSample = FALSE)
NULL



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



#' @rdname plotDot
#' @export
setMethod(
    "plotDot",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        perSample = TRUE,
        colMin = -2.5,
        colMax = 2.5,
        dotMin = 0L,
        dotScale = 6L,
        color = getOption("pointillism.discrete.color", NULL),
        legend = getOption("pointillism.legend", TRUE),
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
)



#' @rdname plotDot
#' @export
setMethod(
    "plotDot",
    signature("seurat"),
    getMethod("plotDot", "SingleCellExperiment")
)

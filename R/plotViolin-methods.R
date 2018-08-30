#' Plot Violin
#'
#' @name plotViolin
#' @family Plot Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams ggplot2::geom_violin
#'
#' @seealso [Seurat::VlnPlot()].
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
#' plotViolin(object, genes = geneIDs)
#' plotViolin(object, genes = geneNames)
#'
#' # Per sample mode can be disabled.
#' plotViolin(object, genes = geneIDs, perSample = FALSE)
NULL



#' @rdname plotViolin
#' @export
setMethod(
    "plotViolin",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        perSample = TRUE,
        scale = c("count", "width", "area"),
        color = getOption("pointillism.discrete.color", NULL),
        legend = getOption("pointillism.legend", TRUE),
        title = NULL
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
)



#' @rdname plotViolin
#' @export
setMethod(
    "plotViolin",
    signature("seurat"),
    getMethod("plotViolin", "SingleCellExperiment")
)

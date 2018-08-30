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
#' # Plotting with either gene IDs or gene names (symbols) works.
#' geneIDs <- head(rownames(object), n = 4L)
#' geneNames <- head(rowRanges(object)[["geneName"]])
#' glimpse(genes)
#' plotViolin(object, genes = genes)
NULL



#' @rdname plotViolin
#' @export
setMethod(
    "plotViolin",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        interestingGroups,
        scale = c("count", "width", "area"),
        color = getOption("pointillism.discrete.color", NULL),
        legend = getOption("pointillism.legend", TRUE),
        title = NULL
    ) {
        validObject(object)
        assert_is_character(genes)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        if (is.character(interestingGroups)) {
            interestingGroups(object) <- interestingGroups
        }
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
        assert_is_all_of(data, "DataFrame")
        assert_is_subset(
            x = c("geneName", "ident", "sampleName"),
            y = colnames(data)
        )

        # Do we need to visualize multiple samples? (logical)
        multiSample <- unique(length(data[["sampleName"]])) > 1L

        if (isTRUE(multiSample)) {
            x <- "sampleName"
        } else {
            x <- "ident"
        }

        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym(x),
                y = !!sym("logcounts"),
                color = !!sym("interestingGroups")
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
                color = paste(interestingGroups, collapse = ":\n")
            ) +
            facet_grid(
                rows = vars(!!sym("ident")),
                cols = vars(!!sym("geneName")),
                scales = "free_y"
            )

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

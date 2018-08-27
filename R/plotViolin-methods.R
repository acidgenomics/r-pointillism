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
#' genes <- head(rownames(object), n = 4L)
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
        fill = getOption("pointillism.discrete.fill", NULL),
        legend = getOption("pointillism.legend", TRUE)
    ) {
        validObject(object)
        assert_is_subset(genes, rownames(object))
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        scale <- match.arg(scale)
        assert_is_any_of(fill, c("ScaleDiscrete", "character", "NULL"))
        if (is.character(fill)) {
            assert_is_a_string(fill)
        }

        # Fetch the gene expression data.
        data <- .fetchGeneData(
            object = object,
            genes = genes,
            assay = "logcounts",
            gene2symbol = TRUE,
            interestingGroups = interestingGroups
        )
        assert_is_subset(
            x = c("gene", "ident", "sampleName"),
            y = colnames(data)
        )
        multiSample <- unique(length(data[["sampleName"]])) > 1L

        # Ensure genes match the data return
        if (isTRUE(.useGene2symbol(object))) {
            g2s <- gene2symbol(object)
            if (length(g2s)) {
                g2s <- g2s[genes, , drop = FALSE]
                genes <- make.unique(g2s[["geneName"]])
                stopifnot(all(genes %in% data[["gene"]]))
            }
        }

        if (isTRUE(multiSample)) {
            x <- "sampleName"
        } else {
            x <- "ident"
        }

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym(x),
                y = !!sym("logcounts"),
                fill = !!sym("interestingGroups")
            )
        ) +
            geom_violin(
                # Never include a color border.
                color = "black",
                scale = scale,
                adjust = 1L,
                show.legend = legend,
                trim = TRUE
            ) +
            # Note that `scales = free_y` will hide the x-axis for some plots.
            labs(
                title = toString(genes),
                fill = paste(interestingGroups, collapse = ":\n")
            ) +
            facet_grid(
                rows = vars(!!sym("ident")),
                cols = vars(!!sym("gene")),
                scales = "free_y"
            )

        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
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

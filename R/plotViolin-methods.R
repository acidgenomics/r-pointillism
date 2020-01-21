#' @rdname plotCounts
#' @name plotViolin
#' @importFrom acidgenerics plotViolin
#' @usage plotViolin(object, ...)
#' @export
NULL



## Updated 2019-09-03.
`plotViolin,SingleCellExperiment` <-  # nolint
    function(
        object,
        genes,
        assay = c("logcounts", "normcounts"),
        perSample = TRUE,
        scale = c("count", "width", "area"),
        color,
        legend,
        title = NULL
    ) {
        validObject(object)
        assert(
            isCharacter(genes),
            isFlag(perSample),
            isGGScale(color, scale = "discrete", aes = "color", nullOK = TRUE),
            isFlag(legend),
            isString(title, nullOK = TRUE)
        )
        assay <- match.arg(assay)
        scale <- match.arg(scale)
        ## Fetch the gene expression data.
        data <- .fetchGeneData(
            object = object,
            genes = genes,
            assay = assay,
            metadata = TRUE
        )
        ## Handling step for multiple samples, if desired.
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
            colorLabs <- "cluster"
        }
        ## Plot.
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym(x),
                y = !!sym(assay),
                color = !!sym(colorMapping)
            )
        ) +
            geom_jitter(show.legend = legend) +
            geom_violin(
                fill = NA,
                scale = scale,
                adjust = 1L,
                show.legend = legend,
                trim = TRUE
            ) +
            ## Note that `scales = free_y` will hide the x-axis for some plots.
            labs(
                x = NULL,
                color = colorLabs,
                title = title
            )
        ## Handling step for multiple samples, if desired.
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
        ## Color.
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
        ## Return.
        p
    }

formals(`plotViolin,SingleCellExperiment`)[
    c("color", "legend")] <-
    list(
        color = discreteColor,
        legend = legend
    )



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotViolin",
    signature = signature("SingleCellExperiment"),
    definition = `plotViolin,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotViolin,Seurat` <-  # nolint
    `plotViolin,SingleCellExperiment`



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotViolin",
    signature = signature("Seurat"),
    definition = `plotViolin,Seurat`
)



## Updated 2019-08-02.
`plotViolin,cell_data_set` <-  # nolint
    `plotViolin,SingleCellExperiment`



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotViolin",
    signature = signature("cell_data_set"),
    definition = `plotViolin,cell_data_set`
)

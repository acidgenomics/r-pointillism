#' @name plotCounts
#' @aliases plotDot plotViolin
#' @inherit bioverbs::plotCounts
#' @note Updated 2019-07-31.
#'
#' @description Visualize genes on a dot or violin plot.
#'
#' @inheritParams acidplots::params
#' @inheritParams basejump::params
#' @inheritParams ggplot2::geom_violin
#' @param colMin `numeric(1)`.
#'   Minimum scaled average expression threshold. Everything smaller will be
#'   set to this.
#' @param colMax `numeric(1)`.
#'   Maximum scaled average expression threshold. Everything larger will be set
#'   to this.
#' @param dotMin `numeric(1)`.
#'   The fraction of cells at which to draw the smallest dot. All cell groups
#'   with less than this expressing the given gene will have no dot drawn.
#' @param dotScale `numeric(1)`.
#'   Scale the size of the points, similar to `cex`.
#' @param geom `character(1)`.
#'   Plot type. Uses [`match.arg()`][base::match.arg] to pick the type.
#'   Currently supports `"dot"` and `"violin"`.
#' @param ... Additional arguments.
#'
#' @seealso
#' - [Seurat::DotPlot()].
#' - [Seurat::VlnPlot()].
#' - [Seurat::RidgePlot()].
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#'
#' ## Plotting with either gene IDs or gene names (symbols) works.
#' genes <- head(rownames(object), n = 4L)
#' print(genes)
#'
#' ## Per sample mode enabled.
#' plotCounts(object, genes = genes, perSample = TRUE)
#'
#' ## Per sample mode disabled.
#' plotCounts(object, genes = genes, perSample = FALSE)
NULL



#' @rdname plotCounts
#' @name plotCounts
#' @importFrom bioverbs plotCounts
#' @usage plotCounts(object, ...)
#' @export
NULL

#' @rdname plotCounts
#' @name plotDot
#' @importFrom bioverbs plotDot
#' @usage plotDot(object, ...)
#' @export
NULL

#' @rdname plotCounts
#' @name plotViolin
#' @importFrom bioverbs plotViolin
#' @usage plotViolin(object, ...)
#' @export
NULL



## plotCounts ==================================================================
## Updated 2019-07-31.
`plotCounts,SingleCellExperiment` <-  # nolint
    function(
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

formals(`plotCounts,SingleCellExperiment`)[["legend"]] <- legend



## plotDot =====================================================================
#' Min Max
#' @note Updated 2019-07-31.
#' @seealso `Seurat:::MinMax`.
#' @noRd
.minMax <- function(data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    data2
}



#' Percent Above
#' @note Updated 2019-07-31.
#' @seealso `Seurat:::PercentAbove`.
#' @noRd
.percentAbove <- function(x, threshold) {
    length(x[x > threshold]) / length(x)
}



## Updated 2019-07-31.
`plotDot,SingleCellExperiment` <-  # nolint
    function(
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
        assert(
            .hasClusters(object),
            isCharacter(genes),
            isFlag(perSample),
            isNumber(colMin),
            isNumber(colMax),
            isNumber(dotMin),
            isNumber(dotScale),
            isGGScale(
                x = color,
                scale = "continuous",
                aes = "colour",
                nullOK = TRUE
            ),
            isFlag(legend),
            isString(title, nullOK = TRUE)
        )

        ## Fetch the gene expression data.
        data <- .fetchGeneData(
            object = object,
            genes = genes,
            assay = "logcounts",
            metadata = TRUE
        )

        ## Prepare data for ggplot.
        ## FIXME Rework metrics return to include ident column for monocle3.
        cols <- c("geneName", "sampleName", "ident")
        data <- data %>%
            as_tibble() %>%
            group_by(!!!syms(cols)) %>%
            summarize(
                avgExp = mean(expm1(!!sym("logcounts"))),
                ## Consider making threshold user definable.
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

        ## Apply our `dotMin` threshold.
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

        ## Handling step for multiple samples, if desired.
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

formals(`plotDot,SingleCellExperiment`)[c(
    "color",
    "legend"
)] <- list(
    color = continuousColorPurpleOrange,
    legend = legend
)



## Updated 2019-07-31.
`plotViolin,SingleCellExperiment` <-  # nolint
    function(
        object,
        genes,
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
            isGGScale(color, scale = "discrete", aes = "colour", nullOK = TRUE),
            isFlag(legend),
            isString(title, nullOK = TRUE)
        )
        scale <- match.arg(scale)

        ## Fetch the gene expression data.
        data <- .fetchGeneData(
            object = object,
            genes = genes,
            assay = "logcounts",
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
            geom_jitter(show.legend = legend) +
            geom_violin(
                fill = NA,
                scale = scale,
                adjust = 1L,
                show.legend = legend,
                trim = TRUE
            ) +
            ## Note that `scales = free_y` will hide the x-axis for some plots.
            labs(title = title, color = colorLabs)

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

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        p
    }

formals(`plotViolin,SingleCellExperiment`)[c("color", "legend")] <-
    list(color = discreteColor, legend = legend)



## Methods =====================================================================
#' @rdname plotCounts
#' @export
setMethod(
    f = "plotCounts",
    signature = signature("SingleCellExperiment"),
    definition = `plotCounts,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotCounts,Seurat` <-  # nolint
    `plotCounts,SingleCellExperiment`



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotCounts",
    signature = signature("Seurat"),
    definition = `plotCounts,SingleCellExperiment`
)



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotDot",
    signature = signature("SingleCellExperiment"),
    definition = `plotDot,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotDot,Seurat` <-  # nolint
    `plotDot,SingleCellExperiment`



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotDot",
    signature = signature("Seurat"),
    definition = `plotDot,Seurat`
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

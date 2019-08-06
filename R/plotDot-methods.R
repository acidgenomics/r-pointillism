#' @rdname plotCounts
#' @name plotDot
#' @importFrom bioverbs plotDot
#' @usage plotDot(object, ...)
#' @export
NULL



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



## Updated 2019-08-06.
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
        value <- "logcounts"

        ## Fetch the gene expression data.
        data <- .fetchGeneData(
            object = object,
            genes = genes,
            value = value,
            metadata = TRUE
        )

        ## Prepare data for ggplot.
        cols <- c("geneName", "sampleName", "ident")
        data <- data %>%
            as_tibble() %>%
            group_by(!!!syms(cols)) %>%
            summarize(
                ## Assuming use of logcounts here.
                avgExp = mean(expm1(!!sym(value))),
                ## Consider making threshold user definable.
                pctExp = .percentAbove(!!sym(value), threshold = 0L)
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
                title = title,
                subtitle = value,
                x = "gene",
                y = "cluster"
            )

        ## Handling step for multiple samples, if desired.
        if (
            isTRUE(perSample) &&
            isTRUE(.hasMultipleSamples(object))
        ) {
            p <- p + facet_wrap(facets = vars(!!sym("sampleName")))
        }

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        p
    }

formals(`plotDot,SingleCellExperiment`)[
    c("color","legend")] <- list(
        color = continuousColorPurpleOrange,
        legend = legend
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



## Updated 2019-08-32.
`plotDot,cell_data_set` <-  # nolint
    `plotDot,SingleCellExperiment`



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotDot",
    signature = signature("cell_data_set"),
    definition = `plotDot,cell_data_set`
)

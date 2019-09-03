#' @rdname plotCounts
#' @name plotDot
#' @importFrom bioverbs plotDot
#' @usage plotDot(object, ...)
#' @export
NULL



#' Min max
#' @note Updated 2019-08-03.
#' @seealso `Seurat:::MinMax`.
#' @noRd
.minMax <- function(x, min, max) {
    x[x > max] <- max
    x[x < min] <- min
    x
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
        x <- .fetchGeneData(
            object = object,
            genes = genes,
            value = value,
            metadata = TRUE
        )
        cols <- c("geneName", "sampleName", "ident")
        assert(isSubset(cols, colnames(x)))
        f <- .group(x[, cols])
        x <- split(x = x, f = f)
        x <- SplitDataFrameList(lapply(
            X = x,
            FUN = function(x) {
                value <- x[["value"]]
                ## Assuming use of logcounts here.
                avgExp <- mean(expm1(value))
                ## Consider making the threshold user definable.
                pctExp <- .percentAbove(value, threshold = 0L)
                DataFrame(
                    geneName = x[["geneName"]][[1L]],
                    sampleName = x[["sampleName"]][[1L]],
                    ident = x[["ident"]][[1L]],
                    avgExp = avgExp,
                    pctExp = pctExp
                )
            }
        ))
        x <- unlist(x, recursive = FALSE, use.names = FALSE)
        ## Calculate the average expression scale per gene.
        x <- split(x, f = x[["geneName"]])
        x <- SplitDataFrameList(lapply(
            X = x,
            FUN = function(x) {
                avgExpScale <- scale(x[["avgExp"]])
                avgExpScale <- .minMax(
                    x = avgExpScale,
                    max = colMax,
                    min = colMin
                )
                x[["avgExpScale"]] <- as.numeric(avgExpScale)
                x
            }
        ))
        x <- unlist(x, recursive = FALSE, use.names = FALSE)
        ## Apply our `dotMin` threshold.
        x[["pctExp"]][x[["pctExp"]] < dotMin] <- NA
        ## Plot.
        data <- as_tibble(x, rownames = NULL)
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
        ## Color.
        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }
        ## Return.
        p
    }

formals(`plotDot,SingleCellExperiment`)[
    c("color", "legend")] <- list(
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

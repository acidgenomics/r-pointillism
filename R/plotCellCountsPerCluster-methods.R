#' @name plotCellCountsPerCluster
#' @inherit AcidGenerics::plotCellCountsPerCluster
#' @note Updated 2020-01-30.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @return Show graphical output. Invisibly return `ggplot`.
#'
#' @examples
#' data(Seurat, package = "AcidTest")
#'
#' ## Seurat ====
#' object <- Seurat
#' plotCellCountsPerCluster(object)
NULL



## Updated 2019-08-03.
`plotCellCountsPerCluster,SingleCellExperiment` <-  # nolint
    function(
        object,
        interestingGroups = NULL
    ) {
        validObject(object)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- cellCountsPerCluster(object = object)
        if (isTRUE(length(levels(data[["sampleName"]])) > 1L)) {
            multipleSamples <- TRUE
            col <- "interestingGroups"
            legendTitle <- paste(interestingGroups, collapse = ":\n")
            showLegend <- TRUE
            xLab <- NULL
        } else {
            multipleSamples <- FALSE
            col <- "ident"
            legendTitle <- NA
            showLegend <- FALSE
            xLab <- "cluster"
        }
        ## Plot.
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym(col),
                y = !!sym("n"),
                fill = !!sym(col)
            )
        ) +
            geom_bar(
                stat = "identity",
                show.legend = showLegend
            ) +
            labs(
                x = xLab,
                y = "n cells",
                fill = legendTitle
            )
        ## Wrap for multiple samples.
        if (isTRUE(multipleSamples)) {
            p <- p + facet_wrap(facets = sym("ident"))
        }
        ## Return.
        p
    }



#' @rdname plotCellCountsPerCluster
#' @export
setMethod(
    f = "plotCellCountsPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = `plotCellCountsPerCluster,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotCellCountsPerCluster,Seurat` <-  # nolint
    `plotCellCountsPerCluster,SingleCellExperiment`



#' @rdname plotCellCountsPerCluster
#' @export
setMethod(
    f = "plotCellCountsPerCluster",
    signature = signature("Seurat"),
    definition = `plotCellCountsPerCluster,Seurat`
)

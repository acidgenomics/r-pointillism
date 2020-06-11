#' @name plotStackedBarPlot
#' @inherit acidgenerics::plotStackedBarPlot
#' @note Updated 2020-06-10.
#'
#' @inheritParams acidroxygen::params
#' @param absolute `logical(1)`.
#'   Return absolute (`TRUE`) or relative/proportional (`FALSE`) cell count.
#' @param ... Additional arguments.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' plotStackedBarPlot(object)
NULL



#' @rdname plotStackedBarPlot
#' @name plotStackedBarPlot
#' @importFrom acidgenerics plotStackedBarPlot
#' @usage plotStackedBarPlot(object, ...)
#' @export
NULL



`plotStackedBarPlot,SingleCellExperiment` <-  # nolint
    function(
        object,
        absolute = FALSE,
        interestingGroups = NULL,
        labels = NULL
    ) {
        validObject(object)
        assert(isFlag(absolute))
        labels <- matchLabels(
            labels = labels,
            choices = eval(formals()[["labels"]])
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        data <- group_by(data, !!!syms(c("interestingGroups", "ident")))
        data <- summarize(data, n = n(), .groups = "keep")
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("interestingGroups"),
                y = !!sym("n"),
                fill = !!sym("ident")
            )
        ) +
            geom_bar(
                color = "black",
                position = ifelse(
                    test = isTRUE(absolute),
                    yes = "stack",
                    no = "fill"
                ),
                stat = "identity"
            )
        if (!isSubset("x", names(labels))) {
            labels[["x"]] <- paste(interestingGroups, collapse = ":\n")
        }
        if (!isSubset("y", names(labels))) {
            labels[["y"]] <- "cell count"
            if (!isTRUE(absolute)) {
                labels[["y"]] <- paste("relative", labels[["y"]])
            }
        }
        p <- p + do.call(what = labs, args = labels)
        p
    }



#' @rdname plotStackedBarPlot
#' @export
setMethod(
    f = "plotStackedBarPlot",
    signature = signature("SingleCellExperiment"),
    definition = `plotStackedBarPlot,SingleCellExperiment`
)



`plotStackedBarPlot,Seurat` <-  # nolint
    `plotStackedBarPlot,SingleCellExperiment`



#' @rdname plotStackedBarPlot
#' @export
setMethod(
    f = "plotStackedBarPlot",
    signature = signature("Seurat"),
    definition = `plotStackedBarPlot,Seurat`
)

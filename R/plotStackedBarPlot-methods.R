#' @name plotStackedBarPlot
#' @inherit AcidGenerics::plotStackedBarPlot
#' @note Updated 2021-03-04.
#'
#' @inheritParams AcidRoxygen::params
#' @param absolute `logical(1)`.
#'   Return absolute (`TRUE`) or relative/proportional (`FALSE`) cell count.
#' @param ... Additional arguments.
#'
#' @examples
#' data(Seurat, package = "AcidTest")
#'
#' ## Seurat ====
#' object <- Seurat
#' plotStackedBarPlot(object)
NULL



## Updated 2021-03-04.
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
        data <- metrics(object, return = "DataFrame")
        ## Generate the summary count table to pass to ggplot.
        cols <- c("interestingGroups", "ident")

        data[, cols]

        f <- as.factor(paste(
            data[["interestingGroups"]],
            data[["ident"]],
            sep = ":"
        ))
        split <- split(x = data, f = f)
        assert(is(split, "SplitDataFrameList"))

        strsplit(x = levels(f), split = ":", fixed = TRUE)

        nrow(split)
        ## unknown:0 unknown:1 unknown:2
        ##        36        25        19

        ## Legacy dplyr alternative approach:
        ## > data <- group_by(data, !!!syms(cols))
        ## > data <- summarize(data, n = n(), .groups = "keep")

        # Groups:   interestingGroups, ident [3]
        ## interestingGroups ident     n
        ## <fct>             <fct> <int>
        ## 1 unknown           0        36
        ## 2 unknown           1        25
        ## 3 unknown           2        19

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

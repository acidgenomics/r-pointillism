## FIXME Consider moving this to AcidPlots.



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
`plotStackedBarPlot,SCE` <-  # nolint
    function(
        object,
        absolute = FALSE,
        interestingGroups = NULL,
        labels = NULL
    ) {
        validObject(object)
        assert(isFlag(absolute))
        labels <- matchLabels(labels)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object, return = "DataFrame")
        ## Generate the summary count table to pass to ggplot.
        cols <- c("interestingGroups", "ident")
        data <- data[, cols, drop = FALSE]
        ## See also our `uniteInterestingGroups()` method, which uses a
        ## similar approach internally.
        f <- as.factor(apply(
            X = as.data.frame(data),
            MARGIN = 1L,
            FUN = paste,
            collapse = ":"
        ))
        split <- split(x = data, f = f)
        assert(is(split, "SplitDataFrameList"))
        data <- as.data.frame(do.call(
            what = rbind,
            args = strsplit(x = levels(f), split = ":", fixed = TRUE)
        ))
        colnames(data) <- cols
        data[["n"]] <- unname(nrow(split))
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
    signature = signature(object = "SingleCellExperiment"),
    definition = `plotStackedBarPlot,SCE`
)

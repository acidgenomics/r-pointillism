#' @name plotPCElbow
#' @inherit AcidGenerics::plotPCElbow
#' @note Updated 2022-06-09.
#'
#' @details
#' Automatically return the smallest number of PCs that match the `minSD`,
#' `minPct`, and `maxCumPct` cutoffs.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @param minPct `numeric(1)` (`0`-`1`).
#' Minimum percent standard deviation.
#'
#' @param maxCumPct `numeric(1)` (`0`-`1`).
#' Maximum cumulative percen standard deviation.
#'
#' @return `ggplot`.
#' Elbow point is saved in `attr(object, "elbow")`.
#'
#' @seealso
#' - `Seurat::PCElbowPlot()`.
#' - `monocle3::plot_pc_variance_explained()`.
#'
#' @examples
#' data(Seurat, package = "AcidTest")
#'
#' ## Seurat ====
#' object <- Seurat
#' plotPCElbow(object)
NULL



## Updated 2019-09-03.
.plotPCElbow <-
    function(pctStdDev,
             minPct = 0.01,
             maxCumPct = 0.9) {
        assert(
            is.numeric(pctStdDev),
            isNumber(minPct),
            allAreInLeftOpenRange(
                x = c(minPct, maxCumPct),
                lower = 0L,
                upper = 1L
            )
        )
        cumsum <- cumsum(pctStdDev)
        xLab <- "PCA component"
        data <- DataFrame(
            pc = seq_along(pctStdDev),
            pct = pctStdDev,
            cumsum = cumsum
        )
        minPctCutoff <- max(data[data[["pct"]] >= minPct, "pc"])
        maxCumPctCutoff <- max(data[data[["cumsum"]] <= maxCumPct, "pc"])
        ## Pick the smallest value of the cutoffs.
        cutoff <- min(minPctCutoff, maxCumPctCutoff)
        ## Percent standard deviation.
        ggpct <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym("pc"),
                y = !!sym("pct")
            )
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = minPct
            ) +
            geom_line() +
            geom_point() +
            geom_vline(xintercept = cutoff) +
            labs(
                x = xLab,
                y = "% std dev"
            ) +
            expand_limits(y = 0L) +
            scale_y_continuous(labels = percent)
        ## Cumulative percent standard deviation.
        ggcumsum <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym("pc"),
                y = !!sym("cumsum")
            )
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = maxCumPct
            ) +
            geom_line() +
            geom_point() +
            geom_vline(xintercept = cutoff) +
            labs(
                x = xLab,
                y = "cum % std dev"
            ) +
            expand_limits(y = c(0L, 1L)) +
            scale_y_continuous(labels = percent)
        ## Return.
        plotlist <- list(
            pct = ggpct,
            cumsum = ggcumsum
        )
        p <- wrap_plots(plotlist)
        attr(p, "elbow") <- cutoff
        p
    }



## Updated 2019-08-02.
`plotPCElbow,Seurat` <- # nolint
    function(object, minPct, maxCumPct) {
        stdDev <- Stdev(object = object, reduction = "pca")
        assert(is.numeric(stdDev))
        pctStdDev <- stdDev^2L / sum(stdDev^2L)
        .plotPCElbow(
            pctStdDev = pctStdDev,
            minPct = minPct,
            maxCumPct = maxCumPct
        )
    }

args <- c("minPct", "maxCumPct")
formals(`plotPCElbow,Seurat`)[args] <- # nolint
    formals(.plotPCElbow)[args]
rm(args)



## > ## Updated 2019-08-02.
## > `plotPCElbow,cell_data_set` <-  # nolint
## >     function(object, minPct, maxCumPct) {
## >         pctStdDev <- slot(object, "preprocess_aux")[["prop_var_expl"]]
## >         .plotPCElbow(
## >             pctStdDev = pctStdDev,
## >             minPct = minPct,
## >             maxCumPct = maxCumPct
## >         )
## >     }
## >
## > args <- c("minPct", "maxCumPct")
## > formals(`plotPCElbow,cell_data_set`)[args] <- formals(.plotPCElbow)[args]



#' @rdname plotPCElbow
#' @export
setMethod(
    f = "plotPCElbow",
    signature = signature(object = "Seurat"),
    definition = `plotPCElbow,Seurat`
)

## > #' @rdname plotPCElbow
## > #' @export
## > setMethod(
## >     f = "plotPCElbow",
## >     signature = signature(object = "cell_data_set"),
## >     definition = `plotPCElbow,cell_data_set`
## > )

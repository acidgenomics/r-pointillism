## FIXME Need to add monocle3 support.



#' @name plotPCElbow
#' @inherit bioverbs::plotPCElbow
#' @note Updated 2019-08-02.
#'
#' @details
#' Automatically return the smallest number of PCs that match the `minSD`,
#' `minPct`, and `maxCumPct` cutoffs.
#'
#' @inheritParams acidplots::params
#' @inheritParams basejump::params
#' @inheritParams params
#' @param minSD `numeric(1)`.
#'   Minimum standard deviation.
#' @param minPct `numeric(1)` (`0`-`1`).
#'   Minimum percent standard deviation.
#' @param maxCumPct `numeric(1)` (`0`-`1`).
#'   Maximum cumulative percen standard deviation.
#' @param ... Additional arguments.
#'
#' @return
#' - Show graphical output of elbow plots.
#' - Invisibly return numeric sequence vector of PCs to include for
#'   dimensionality reduction analysis.
#'
#' @seealso
#' - `Seurat::PCElbowPlot()`.
#' - `monocle3::plot_pc_variance_explained()`.
#'
#' @examples
#' data(
#'     Seurat,
#'     cell_data_set,
#'     package = "acidtest"
#' )
#'
#' ## Seurat ====
#' object <- Seurat
#' plotPCElbow(object)
#'
#' ## cell_data_set ====
#' object <- cell_data_set
#' plotPCElbow(object)
NULL



#' @rdname plotPCElbow
#' @name plotPCElbow
#' @importFrom bioverbs plotPCElbow
#' @usage plotPCElbow(object, ...)
#' @export
NULL



## Updated 2019-07-31.
`plotPCElbow,Seurat` <-  # nolint
    function(
        object,
        reduction = "PCA",
        minSD = 1L,
        minPct = 0.01,
        maxCumPct = 0.9,
        trans = c("identity", "sqrt")
    ) {
        assert(
            isString(reduction),
            isNumber(minSD),
            isPositive(minSD),
            isNumber(minPct),
            isNumber(maxCumPct),
            allAreInLeftOpenRange(
                x = c(minPct, maxCumPct),
                lower = 0L,
                upper = 1L
            )
        )
        trans <- match.arg(trans)
        xLab <- "PCA component"

        ## dr: dimensional reduction
        ## sdev: standard deviation
        ## FIXME Extract the code for this...
        ## FIXME Seurat:::Stdev.DimReduc
        ## https://github.com/satijalab/seurat/blob/master/R/objects.R#L3386
        sdev <- Stdev(object = object, reduction = tolower(reduction))
        assert(is.numeric(sdev))
        pct <- sdev ^ 2L / sum(sdev ^ 2L)
        cumsum <- cumsum(pct)

        data <- tibble(
            pc = seq_along(sdev),
            sdev = sdev,
            pct = pct,
            cumsum = cumsum
        )

        minSDCutoff <- data %>%
            .[.[["sdev"]] >= minSD, "pc"] %>%
            max()
        minPctCutoff <- data %>%
            .[.[["pct"]] >= minPct, "pc"] %>%
            max()
        maxCumPctCutoff <- data %>%
            .[.[["cumsum"]] <= maxCumPct, "pc"] %>%
            max()

        ## Pick the smallest value of the cutoffs
        cutoff <- min(minSDCutoff, minPctCutoff, maxCumPctCutoff)

        ## Standard deviation --------------------------------------------------
        ggsd <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("pc"),
                y = !!sym("sdev")
            )
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = minSD
            ) +
            geom_line() +
            geom_point() +
            geom_vline(xintercept = cutoff) +
            labs(
                x = xLab,
                y = "std dev"
            ) +
            expand_limits(y = 0L) +
            scale_y_continuous(trans = trans)

        ## Percent standard deviation ------------------------------------------
        ggpct <- ggplot(
            data = data,
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
            scale_y_continuous(labels = percent, trans = trans)

        ## Cumulative percent standard deviation -------------------------------
        ggcumsum <- ggplot(
            data = data,
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
            scale_y_continuous(labels = percent, trans = trans)

        plotlist <- list(
            sd = ggsd,
            pct = ggpct,
            cumsum = ggcumsum
        )

        p <- plot_grid(plotlist = plotlist)
        show(p)

        invisible(seq_len(cutoff))
    }



#' @rdname plotPCElbow
#' @export
setMethod(
    f = "plotPCElbow",
    signature = signature("Seurat"),
    definition = `plotPCElbow,Seurat`
)



`plotPCElbow,cell_data_set` <-  # nolint
    function(object) {
        ## FIXME Rework this to be consistent with Seurat method.
        ## FIXME Probably need to break out the plot code.
        monocle3::plot_pc_variance_explained(object)
    }



#' @rdname plotPCElbow
#' @export
setMethod(
    f = "plotPCElbow",
    signature = signature("cell_data_set"),
    definition = `plotPCElbow,cell_data_set`
)

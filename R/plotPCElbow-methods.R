#' @name plotPCElbow
#' @inherit bioverbs::plotPCElbow
#' @inheritParams basejump::params
#'
#' @details
#' Automatically return the smallest number of PCs that match the `minSD`,
#' `minPct`, and `maxCumPct` cutoffs.
#'
#' @param minSD `numeric(1)`.
#'   Minimum standard deviation.
#' @param minPct `numeric(1)` (`0`-`1`).
#'   Minimum percent standard deviation.
#' @param maxCumPct `numeric(1)` (`0`-`1`).
#'   Maximum cumulative percen standard deviation.
#'
#' @return
#' - Show graphical output of elbow plots.
#' - Invisibly return numeric sequence vector of PCs to include for
#'   dimensionality reduction analysis.
#'
#' @seealso `Seurat::PCElbowPlot()`.
#'
#' @examples
#' data(pbmc_small, package = "Seurat")
#' plotPCElbow(pbmc_small)
NULL



#' @importFrom bioverbs plotPCElbow
#' @aliases NULL
#' @export
bioverbs::plotPCElbow



plotPCElbow.seurat <-  # nolint
    function(
        object,
        minSD = 1L,
        minPct = 0.01,
        maxCumPct = 0.9,
        trans = c("identity", "sqrt")
    ) {
        assert(
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

        # dr: dimensional reduction
        # sdev: standard deviation
        sdev <- object@dr[["pca"]]@sdev
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

        # Pick the smallest value of the cutoffs
        cutoff <- min(minSDCutoff, minPctCutoff, maxCumPctCutoff)

        # Standard deviation ---------------------------------------------------
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
                x = "PC",
                y = "std dev"
            ) +
            expand_limits(y = 0L) +
            scale_y_continuous(trans = trans)

        # Percent standard deviation -------------------------------------------
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
                x = "PC",
                y = "% std dev"
            ) +
            expand_limits(y = 0L) +
            scale_y_continuous(labels = percent, trans = trans)

        # Cumulative percent standard deviation --------------------------------
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
                x = "PC",
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
    signature = signature("seurat"),
    definition = plotPCElbow.seurat
)

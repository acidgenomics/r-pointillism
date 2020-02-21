#' @name plotFeature
#' @inherit acidgenerics::plotFeature
#' @note Updated 2020-02-21.
#'
#' @inheritParams acidroxygen::params
#' @param features `character`. Features to plot (e.g. gene expression, PC
#'   scores, number of genes detected).
#' @param ... Additional arguments.
#'
#' @seealso [Seurat::FeaturePlot()].
#'
#' @return `ggplot` (1 feature) or `list` (multiple features).
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' plotFeature(
#'     object = object,
#'     features = c("nCount_RNA", "nFeature_RNA", "PC_1", "PC_2"),
#'     reduction = "UMAP"
#' )
NULL



#' @rdname plotFeature
#' @name plotFeature
#' @importFrom acidgenerics plotFeature
#' @usage plotFeature(object, ...)
#' @export
NULL



## Updated 2020-02-21.
`plotFeature,SingleCellExperiment` <-  # nolint
    function(
        object,
        features,
        reduction,
        color,
        pointSize,
        pointAlpha,
        pointsAsNumbers,
        label,
        labelSize,
        dark,
        legend
    ) {
        assert(
            isCharacter(features),
            isScalar(reduction),
            isGGScale(
                x = color,
                scale = "continuous",
                aes = "color",
                nullOK = TRUE
            ),
            isNumber(pointSize),
            isNumber(pointAlpha),
            isFlag(pointsAsNumbers),
            isFlag(label),
            isNumber(labelSize),
            isFlag(dark),
            isFlag(legend)
        )
        ## Dark mode setup.
        if (isTRUE(dark)) {
            fill <- "black"
        } else {
            fill <- "white"
        }
        ## Get the dimension reduction data.
        data <- .fetchReductionData(
            object = object,
            reduction = reduction
        )
        ## Label the axes.
        axes <- colnames(data)[seq_len(2L)]
        ## If the features are not defined, attempt to merge all reduced dims
        ## information before stopping.
        if (!isSubset(features, colnames(data))) {
            reductionData <- .bindReducedDims(object)
            data <- data[, setdiff(colnames(data), colnames(reductionData))]
            data <- cbind(data, reductionData)
        }
        ## Only supporting visualization of numeric column metadata.
        supported <- bapply(X = data, FUN = is.numeric)
        supported <- names(supported)[supported]
        ## These values are required for dim reduction labeling.
        blacklist <- c("centerX", "centerY", "x", "y")
        supported <- setdiff(supported, blacklist)
        if (!isSubset(features, supported)) {
            setdiff <- setdiff(features, supported)
            stop(sprintf(
                fmt = paste0(
                    "%s ",
                    ngettext(
                        n = length(setdiff),
                        msg1 = "feature",
                        msg2 = "features"
                    ),
                    " not defined: %s\n",
                    "Available:\n%s"
                ),
                length(setdiff),
                toString(setdiff, width = 200L),
                printString(supported)
            ))
        }
        plotlist <- lapply(
            X = features,
            FUN = function(feature) {
                p <- ggplot(
                    data = as.data.frame(data),
                    mapping = aes(
                        x = !!sym("x"),
                        y = !!sym("y"),
                        color = !!sym(feature)
                    )
                )
                if (isTRUE(pointsAsNumbers)) {
                    if (pointSize < 4L) pointSize <- 4L
                    p <- p +
                        geom_text(
                            mapping = aes(
                                x = !!sym("x"),
                                y = !!sym("y"),
                                label = !!sym("ident"),
                                color = !!sym(feature)
                            ),
                            alpha = pointAlpha,
                            size = pointSize
                        )
                } else {
                    p <- p +
                        geom_point(
                            alpha = pointAlpha,
                            size = pointSize
                        )
                }
                p <- p +
                    labs(
                        color = NULL,
                        x = axes[[1L]],
                        y = axes[[2L]],
                        title = feature
                    )
                if (isTRUE(label)) {
                    if (isTRUE(dark)) {
                        labelColor <- "white"
                    } else {
                        labelColor <- "black"
                    }
                    p <- p +
                        geom_text(
                            mapping = aes(
                                x = !!sym("centerX"),
                                y = !!sym("centerY"),
                                label = !!sym("ident")
                            ),
                            color = labelColor,
                            size = labelSize,
                            fontface = "bold"
                        )
                }
                if (isTRUE(dark)) {
                    p <- p + acid_theme_dark()
                    if (is.null(color)) {
                        color <- darkMarkerColors
                    }
                }
                if (is(color, "ScaleContinuous")) {
                    p <- p + color
                }
                if (!isTRUE(legend)) {
                    p <- p + guides(color = "none")
                }
                p
            }
        )

        ## Return --------------------------------------------------------------
        if (isTRUE(length(features) > 1L)) {
            plot_grid(plotlist = plotlist) +
                theme(
                    plot.background = element_rect(color = NA, fill = fill)
                )
        } else {
            plotlist[[1L]]
        }
    }

args <- c(
    "dark",
    "expression",
    "label",
    "labelSize",
    "legend",
    "pointAlpha",
    "pointSize",
    "pointsAsNumbers",
    "reduction"
)
args1 <- c(args, "color")
args2 <- c(args, "continuousColor")
formals(`plotFeature,SingleCellExperiment`)[args1] <- .formalsList[args2]
rm(args, args1, args2)



#' @rdname plotFeature
#' @export
setMethod(
    f = "plotFeature",
    signature = signature("SingleCellExperiment"),
    definition = `plotFeature,SingleCellExperiment`
)



## Updated 2019-07-31.
`plotFeature,Seurat` <-  # nolint
    `plotFeature,SingleCellExperiment`



#' @rdname plotFeature
#' @export
setMethod(
    f = "plotFeature",
    signature = signature("Seurat"),
    definition = `plotFeature,Seurat`
)

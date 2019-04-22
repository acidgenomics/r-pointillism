#' @name plotCellCountsPerCluster
#' @inherit bioverbs::plotCellCountsPerCluster
#' @inheritParams acidplots::params
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @return Show graphical output. Invisibly return `ggplot`.
#'
#' @examples
#' data(seurat)
#' plotCellCountsPerCluster(seurat)
NULL



#' @rdname plotCellCountsPerCluster
#' @name plotCellCountsPerCluster
#' @importFrom bioverbs plotCellCountsPerCluster
#' @export
NULL



plotCellCountsPerCluster.SingleCellExperiment <-  # nolint
    function(
        object,
        interestingGroups = NULL
    ) {
        validObject(object)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- cellCountsPerCluster(
            object = object,
            interestingGroups = interestingGroups
        )
        assert(is(data, "tbl_df"))
        ggplot(
            data = data,
            mapping = aes(
                x = !!sym("interestingGroups"),
                y = !!sym("n"),
                fill = !!sym("interestingGroups")
            )
        ) +
            geom_bar(stat = "identity") +
            labs(
                x = NULL,
                y = "n cells",
                fill = paste(interestingGroups, collapse = ":\n")
            ) +
            facet_wrap(facets = sym("ident"))
    }



#' @rdname plotCellCountsPerCluster
#' @export
setMethod(
    f = "plotCellCountsPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = plotCellCountsPerCluster.SingleCellExperiment
)



plotCellCountsPerCluster.Seurat <-  # nolint
    plotCellCountsPerCluster.SingleCellExperiment



#' @rdname plotCellCountsPerCluster
#' @export
setMethod(
    f = "plotCellCountsPerCluster",
    signature = signature("Seurat"),
    definition = plotCellCountsPerCluster.Seurat
)

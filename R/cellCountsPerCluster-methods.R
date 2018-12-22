#' @name cellCountsPerCluster
#' @inherit bioverbs::cellCountsPerCluster
#' @inheritParams basejump::params
#'
#' @return `tbl_df`. Grouped by `ident` column and arranged by `n`.
#'
#' @examples
#' data(seurat_small)
#' x <- cellCountsPerCluster(seurat_small)
#' print(x)
NULL



#' @importFrom bioverbs cellCountsPerCluster
#' @aliases NULL
#' @export
bioverbs::cellCountsPerCluster



cellCountsPerCluster.SingleCellExperiment <-  # nolint
    function(object, interestingGroups = NULL) {
        validObject(object)
        assert(.hasIdent(object))
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        data <- metrics(object)
        cols <- unique(c("ident", interestingGroups))
        assert(isSubset(cols, colnames(data)))
        data %>%
            as_tibble() %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!!syms(cols)) %>%
            summarize(n = n()) %>%
            ungroup() %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!sym("ident")) %>%
            mutate(ratio = !!sym("n") / sum(!!sym("n")))
    }



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    f = "cellCountsPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = cellCountsPerCluster.SingleCellExperiment
)



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    f = "cellCountsPerCluster",
    signature = signature("seurat"),
    definition = cellCountsPerCluster.SingleCellExperiment
)

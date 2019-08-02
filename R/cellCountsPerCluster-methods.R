#' @name cellCountsPerCluster
#' @inherit bioverbs::cellCountsPerCluster
#' @note Updated 2019-07-31.
#'
#' @inheritParams basejump::params
#' @param ... Additional arguments.
#'
#' @return `tbl_df`. Grouped by `ident` column and arranged by `n`.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' x <- cellCountsPerCluster(object)
#' print(x)
NULL



#' @rdname cellCountsPerCluster
#' @name cellCountsPerCluster
#' @importFrom bioverbs cellCountsPerCluster
#' @usage cellCountsPerCluster(object, ...)
#' @export
NULL



## FIXME Need to rework this for monocle3.
## Updated 2019-08-02.
`cellCountsPerCluster,SingleCellExperiment` <-  # nolint
    function(object, interestingGroups = NULL) {
        validObject(object)
        assert(.hasClusters(object))
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        ## FIXME Rework metrics return to include ident column for monocle3.
        data <- metrics(object)
        cols <- unique(c("ident", interestingGroups, "interestingGroups"))
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
    definition = `cellCountsPerCluster,SingleCellExperiment`
)



## Updated 2019-07-31.
`cellCountsPerCluster,Seurat` <-  # nolint
    `cellCountsPerCluster,SingleCellExperiment`



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    f = "cellCountsPerCluster",
    signature = signature("Seurat"),
    definition = `cellCountsPerCluster,Seurat`
)

# FIXME Add a corresponding `plot` function.



#' Cell Counts per Cluster
#'
#' @name cellCountsPerCluster
#' @family Cluster Statistics Functions
#'
#' @inheritParams general
#'
#' @return `tbl_df`. Grouped by `ident` column and arranged by `n`.
#'
#' @examples
#' data(seurat_small)
#' x <- cellCountsPerCluster(seurat_small)
#' print(x)
NULL



.cellCountsPerCluster.SCE <-  # nolint
    function(object, interestingGroups = NULL) {
        validObject(object)
        .assertHasIdent(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        data <- metrics(object)
        cols <- unique(c("ident", interestingGroups))
        assert_is_subset(cols, colnames(data))
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
    definition = .cellCountsPerCluster.SCE
)



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    f = "cellCountsPerCluster",
    signature = signature("seurat"),
    definition = getMethod(
        f = "cellCountsPerCluster",
        signature = signature("SingleCellExperiment")
    )
)

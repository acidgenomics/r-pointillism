#' Cell Counts per Cluster
#'
#' @name cellCountsPerCluster
#' @family Cluster Statistics Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `tibble` grouped by "`ident`" column, arranged by abundance.
#'
#' @examples
#' x <- cellCountsPerCluster(sce_small)
#' glimpse(x)
NULL



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    "cellCountsPerCluster",
    signature("SingleCellExperiment"),
    function(object, interestingGroups) {
        validObject(object)
        .assertHasIdent(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        metrics <- metrics(object, interestingGroups = interestingGroups)
        cols <- unique(c("ident", interestingGroups))
        assert_is_subset(cols, colnames(metrics))
        metrics %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!!syms(cols)) %>%
            summarize(n = n()) %>%
            ungroup() %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!sym("ident")) %>%
            mutate(ratio = !!sym("n") / sum(!!sym("n")))
    }
)



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    "cellCountsPerCluster",
    signature("seurat"),
    getMethod("cellCountsPerCluster", "SingleCellExperiment")
)

#' Cluster Cell Counts per Sample
#'
#' @name clusterCellCountsPerSample
#' @family Cluster Statistics Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `grouped_df`. Grouped by `sampleName` column, arranged by abundance.
#'
#' @examples
#' data(seurat_small)
#' x <- clusterCellCountsPerSample(seurat_small)
#' print(x)
NULL



.clusterCellCountsPerSample.SCE <-  # nolint
    function(object) {
        .assertHasIdent(object)
        metrics <- metrics(object)
        cols <- c("sampleName", "ident")
        assert_is_subset(cols, colnames(metrics))
        metrics %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!!syms(cols)) %>%
            summarize(n = n()) %>%
            ungroup() %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!sym("sampleName")) %>%
            mutate(ratio = !!sym("n") / sum(!!sym("n")))
    }



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    f = "clusterCellCountsPerSample",
    signature = signature("SingleCellExperiment"),
    definition = .clusterCellCountsPerSample.SCE
)



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    f = "clusterCellCountsPerSample",
    signature = signature("seurat"),
    definition = getMethod(
        f = "clusterCellCountsPerSample",
        signature = signature("SingleCellExperiment")
    )
)

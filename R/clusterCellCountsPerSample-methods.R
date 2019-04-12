#' @name clusterCellCountsPerSample
#' @inherit bioverbs::clusterCellCountsPerSample
#' @inheritParams basejump::params
#' @examples
#' data(seurat_small)
#' x <- clusterCellCountsPerSample(seurat_small)
#' print(x)
NULL



#' @importFrom bioverbs clusterCellCountsPerSample
#' @aliases NULL
#' @export
bioverbs::clusterCellCountsPerSample



clusterCellCountsPerSample.SingleCellExperiment <-  # nolint
    function(object) {
        assert(.hasIdent(object))
        metrics <- metrics(object)
        cols <- c("sampleName", "ident")
        assert(isSubset(cols, colnames(metrics)))
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
    definition = clusterCellCountsPerSample.SingleCellExperiment
)



clusterCellCountsPerSample.Seurat <-  # nolint
    clusterCellCountsPerSample.SingleCellExperiment



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    f = "clusterCellCountsPerSample",
    signature = signature("Seurat"),
    definition = clusterCellCountsPerSample.Seurat
)

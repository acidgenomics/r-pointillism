#' @name clusterCellCountsPerSample
#' @inherit bioverbs::clusterCellCountsPerSample
#' @note Updated 2019-07-31.
#'
#' @inheritParams basejump::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' x <- clusterCellCountsPerSample(object)
#' print(x)
NULL



#' @rdname clusterCellCountsPerSample
#' @name clusterCellCountsPerSample
#' @importFrom bioverbs clusterCellCountsPerSample
#' @usage clusterCellCountsPerSample(object, ...)
#' @export
NULL



## Updated 2019-08-02.
`clusterCellCountsPerSample,SingleCellExperiment` <-  # nolint
    function(object) {
        assert(.hasClusters(object))
        ## FIXME Rework metrics return to include ident column for monocle3.
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
    definition = `clusterCellCountsPerSample,SingleCellExperiment`
)



## Updated 2019-07-31.
`clusterCellCountsPerSample,Seurat` <-  # nolint
    `clusterCellCountsPerSample,SingleCellExperiment`



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    f = "clusterCellCountsPerSample",
    signature = signature("Seurat"),
    definition = `clusterCellCountsPerSample,Seurat`
)

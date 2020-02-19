#' @name pseudobulk
#' @inherit acidgenerics::pseudobulk
#' @note Updated 2020-01-30.
#'
#' @inheritParams basejump::aggregate
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(SingleCellExperiment, package = "acidtest")
#'
#' ## SingleCellExperiment ====
#' object <- SingleCellExperiment
#' x <- pseudobulk(object)
#' print(x)
NULL



#' @rdname pseudobulk
#' @name pseudobulk
#' @importFrom acidgenerics pseudobulk
#' @usage pseudobulk(object, ...)
#' @export
NULL



## Updated 2020-01-30.
`pseudobulk,SingleCellExperiment` <-  # nolint
    function(
        object,
        fun = c("sum", "mean")
    ) {
        validObject(object)
        fun <- match.arg(fun)
        rse <- as(object, "RangedSummarizedExperiment")
        colData <- colData(rse)
        assert(areDisjointSets("aggregate", colnames(colData)))
        colData[["aggregate"]] <- cell2sample(object)
        if ("sampleID" %in% colnames(colData)) {
            colData[["sampleID"]] <- NULL
        }
        colData(rse) <- colData
        aggregateCols(x = rse, col = "aggregate", fun = fun)
    }



#' @rdname pseudobulk
#' @export
setMethod(
    f = "pseudobulk",
    signature = signature("SingleCellExperiment"),
    definition = `pseudobulk,SingleCellExperiment`
)

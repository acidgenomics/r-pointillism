#' Size factors
#'
#' @inheritParams basejump::params
#'
#' @name sizeFactors
#'
#' @return `numeric`.
#'   Named numeric vector, corresponding to cells.
#'
#' @examples
#' data(SingleCellExperiment, package = "acidtest")
#'
#' ## SingleCellExperiment ====
#' object <- SingleCellExperiment
#' object <- estimateSizeFactors(object)
#' head(sizeFactors(object))
NULL



## FIXME Need to reexport generic



## Updated 2019-08-02.
`sizeFactors,SingleCellExperiment` <-  # nolint
    methodFunction(
        f = "sizeFactors",
        signature = "DESeqDataSet",
        package = "DESeq2"
    )



#' @rdname sizeFactors
#' @export
setMethod(
    f = "sizeFactors",
    signature = signature("SingleCellExperiment"),
    definition = `sizeFactors,SingleCellExperiment`
)



## Updated 2019-08-02.
`sizeFactors,cell_data_set` <-  # nolint
    function(object) {
        monocle3::size_factors(object)
    }



#' @rdname sizeFactors
#' @export
setMethod(
    f = "sizeFactors",
    signature = signature("cell_data_set"),
    definition = `sizeFactors,cell_data_set`
)

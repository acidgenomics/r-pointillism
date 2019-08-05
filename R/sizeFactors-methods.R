## Consider adding support for: `clearSizeFactors()`, `sizeFactorNames()`.



#' Size factors
#'
#' @name sizeFactors
#' @note Updated 2019-08-05.
#'
#' @inheritParams basejump::params
#' @param value Value to be assigned to corresponding components of object.
#' @param ... Additional arguments.
#'
#' @return `numeric`.
#'   Named numeric vector, corresponding to cells.
#'
#' @seealso
#' - `monocle3::size_factors()`.
#'
#' @examples
#' data(SingleCellExperiment, package = "acidtest")
#'
#' ## SingleCellExperiment ====
#' object <- SingleCellExperiment
#' object <- estimateSizeFactors(object)
#' head(sizeFactors(object))
NULL



#' @rdname sizeFactors
#' @name sizeFactors
#' @importFrom BiocGenerics sizeFactors
#' @usage sizeFactors(object, ...)
#' @export
NULL

#' @rdname sizeFactors
#' @name sizeFactors<-
#' @importFrom BiocGenerics sizeFactors<-
#' @usage sizeFactors(object, ...) <- value
#' @export
NULL



## Updated 2019-08-02.
`sizeFactors,SingleCellExperiment` <-  # nolint
    methodFunction(
        f = "sizeFactors",
        signature = "DESeqDataSet",
        package = "DESeq2"
    )



## Updated 2019-08-05.
`sizeFactors,cell_data_set` <-  # nolint
    function(object, type = NULL) {
        assert(is.null(type))
        colData(object)[["Size_Factor"]]
    }



#' @rdname sizeFactors
#' @export
setMethod(
    f = "sizeFactors",
    signature = signature("cell_data_set"),
    definition = `sizeFactors,cell_data_set`
)



## Updated 2019-08-05.
`sizeFactors<-,cell_data_set,numeric` <-  # nolint
    function(object, value) {
        value <- unname(value)
        colData(object)[["Size_Factor"]] <- value
        object
    }



#' @rdname sizeFactors
#' @export
setReplaceMethod(
    f = "sizeFactors",
    signature = signature(
        object = "cell_data_set",
        value = "numeric"
    ),
    definition = `sizeFactors<-,cell_data_set,numeric`
)

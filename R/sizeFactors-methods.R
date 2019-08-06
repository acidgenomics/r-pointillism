## FIXME SingleCellExperiment is currently putting in colData instead of intColData...



#' @name sizeFactors
#' @inherit basejump::sizeFactors
#'
#' @note For `Seurat` objects, use the `Seurat::NormalizeData()` function or the
#'   sctransform package to normalize counts, rather than using size factors.
#' @note Updated 2019-08-06.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @seealso
#' - `monocle3::size_factors()`.
#'
#' @examples
#' data(
#'     SingleCellExperiment,
#'     cell_data_set,
#'     package = "acidtest"
#' )
#'
#' ## SingleCellExperiment ====
#' object <- SingleCellExperiment
#' object <- estimateSizeFactors(object)
#' head(sizeFactors(object))
#'
#' ## cell_data_set ====
#' object <- cell_data_set
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



## Updated 2019-08-06.
`sizeFactors,cell_data_set` <-  # nolint
    function(object) {
        colData(object)[["Size_Factor"]]
    }



#' @rdname sizeFactors
#' @export
setMethod(
    f = "sizeFactors",
    signature = signature("cell_data_set"),
    definition = `sizeFactors,cell_data_set`
)



## nolint start
##
## > getMethod(
## >     f = "sizeFactors<-",
## >     signature = signature(
## >         object = "SummarizedExperiment",
## >         value = "numeric"
## >     ),
## >     where = asNamespace("basejump")
## > )
##
## nolint end



## Updated 2019-08-06.
`sizeFactors<-,cell_data_set,numeric` <-  # nolint
    function(object, value) {
        assert(
            all(!is.na(value)),
            all(is.finite(value)),
            all(value > 0)
        )
        colData(object)[["Size_Factor"]] <- unname(value)
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

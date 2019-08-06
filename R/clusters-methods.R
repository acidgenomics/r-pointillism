#' @name clusters
#' @inherit bioverbs::clusters
#' @note Updated 2019-08-02.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @seealso
#' - `Seurat::Idents()`.
#' - `monocle3::clusters()`.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' x <- clusters(object)
#' head(x)
#' table(x)
NULL



#' @rdname clusters
#' @name clusters
#' @importFrom bioverbs clusters
#' @usage clusters(object, ...)
#' @export
NULL



## Updated 2019-08-06.
`clusters,SingleCellExperiment` <-  # nolint
    function(object) {
        validObject(object)
        assert(isSubset("ident", colnames(colData(object))))
        ident <- colData(object)[["ident"]]
        assert(is.factor(ident))
        names(ident) <- colnames(object)
        ident
    }



#' @rdname clusters
#' @export
setMethod(
    f = "clusters",
    signature = signature("SingleCellExperiment"),
    definition = `clusters,SingleCellExperiment`
)



## Updated 2019-08-02.
`clusters,Seurat` <-  # nolint
    function(object) {
        validObject(object)
        Idents(object)
    }



#' @rdname clusters
#' @export
setMethod(
    f = "clusters",
    signature = signature("Seurat"),
    definition = `clusters,Seurat`
)



## Note that monocle3 currently defines the generic using "x" instead of
## "object", and requires "reduction_method" in definition.
##
## For `reduction`, note that positional scalar works.
##
## Updated 2019-08-02.
`clusters,cell_data_set` <-  # nolint
    function(object, reduction) {
        validObject(object)
        assert(isScalar(reduction))
        monocle3::clusters(
            x = object,
            reduction_method = reduction
        )
    }

f <- methodFormals(
    f = "clusters",
    signature = signature(x = "cell_data_set"),
    package = "monocle3"
)
formals(`clusters,cell_data_set`)[["reduction"]] <-
    f[["reduction_method"]]



#' @rdname clusters
#' @export
setMethod(
    f = "clusters",
    signature = signature("cell_data_set"),
    definition = `clusters,cell_data_set`
)

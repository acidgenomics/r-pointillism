#' @name clusters
#' @inherit bioverbs::clusters
#' @note Updated 2019-08-02.
#'
#' @inheritParams basejump::params
#' @param ... Additional arguments.
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
## For `reductionMethod`, note that positional scalar works.
##
## Updated 2019-08-02.
`clusters,cell_data_set` <-  # nolint
    function(object, reductionMethod) {
        validObject(object)
        assert(isScalar(reductionMethod))

        monocle3::clusters(
            x = object,
            reduction_method = reductionMethod
        )
    }

f <- methodFormals(
    f = "clusters",
    signature = signature(x = "cell_data_set"),
    package = "monocle3"
)
formals(`clusters,cell_data_set`)[["reductionMethod"]] <-
    f[["reduction_method"]]



#' @rdname clusters
#' @export
setMethod(
    f = "clusters",
    signature = signature("cell_data_set"),
    definition = `clusters,cell_data_set`
)



## Consider removing this fallback method support.
## Updated 2019-08-02.
`clusters,SingleCellExperiment` <-  # nolint
    function(object) {
        validObject(object)

        ## Look for `ident` column defined by Seurat.
        if (isSubset("ident", colnames(colData(object)))) {
            x <- colData(object)[["ident"]]
            assert(is.factor(x))
            names(x) <- colnames(object)
            return(x)
        }

        ## Look for `clusters` slot defined by monocle3.
        if (.hasSlot("clusters", object)) {
            x <- slot(object, "clusters")
            ## Here we're assuming the first slot, typically UMAP.
            x <- x[[1L]]
            assert(isSubset("clusters", names(x)))
            x <- x[["clusters"]]
            assert(
                is.factor(x),
                identical(names(x), colnames(object))
            )
            return(x)
        }

        stop("Failed to locate cluster mappings in the object.")
    }



#' @rdname clusters
#' @export
setMethod(
    f = "clusters",
    signature = signature("SingleCellExperiment"),
    definition = `clusters,SingleCellExperiment`
)

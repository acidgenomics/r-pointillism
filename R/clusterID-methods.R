## FIXME Consider renaming to clusters from clusterID.



#' @name clusterID
#' @inherit bioverbs::clusterID
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
#' x <- clusterID(object)
#' head(x)
#' table(x)
NULL



#' @rdname clusterID
#' @name clusterID
#' @importFrom bioverbs clusterID
#' @usage clusterID(object, ...)
#' @export
NULL



## Updated 2019-08-02.
`clusterID,Seurat` <-  # nolint
    function(object) {
        validObject(object)
        Idents(object)
    }



#' @rdname clusterID
#' @export
setMethod(
    f = "clusterID",
    signature = signature("Seurat"),
    definition = `clusterID,Seurat`
)



## For `reductionMethod`, note that positional scalar works.
## Updated 2019-08-02.
`clusterID,cell_data_set` <-  # nolint
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
formals(`clusterID,cell_data_set`)[["reductionMethod"]] <-
    f[["reduction_method"]]



#' @rdname clusterID
#' @export
setMethod(
    f = "clusterID",
    signature = signature("cell_data_set"),
    definition = `clusterID,cell_data_set`
)



## Consider removing this fallback method support.
## Updated 2019-08-02.
`clusterID,SingleCellExperiment` <-  # nolint
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



#' @rdname clusterID
#' @export
setMethod(
    f = "clusterID",
    signature = signature("SingleCellExperiment"),
    definition = `clusterID,SingleCellExperiment`
)

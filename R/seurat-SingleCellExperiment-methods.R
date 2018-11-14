#' Extend S4 Methods for `seurat` Class
#'
#' Provide limited `SingleCellExperiment`-like method support.
#'
#' @name seurat-SingleCellExperiment
#' @keywords internal
#'
#' @inheritParams basejump::params
#'
#' @return Match `SummarizedExperiment` method return.
NULL



# Internal =====================================================================
.getSeuratStash <- function(object, name) {
    assert_that(is(object, "seurat"))
    assert_is_a_string(name)

    misc <- slot(object, "misc")

    # Early return if the `misc` slot is `NULL`.
    if (is.null(misc)) {
        return(NULL)
    }

    # Look first directly in `object@misc` slot.
    x <- misc[[name]]
    if (!is.null(x)) {
        return(x)
    }

    # Next, handle legacy `bcbio` stash list inside `object@misc`.
    # As of v0.1.3, stashing directly into `object@misc`.
    if ("bcbio" %in% names(misc)) {
        x <- misc[["bcbio"]][[name]]
        if (!is.null(x)) {
            return(x)
        }
    }

    NULL
}



# assay ========================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(
    f = "assay",
    signature = signature("seurat"),
    definition = function(x, ...) {
        assay(as.SingleCellExperiment(x), ...)
    }
)



# assayNames ===================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assayNames
#' @export
setMethod(
    f = "assayNames",
    signature = signature("seurat"),
    definition = function(x, ...) {
        assayNames(as.SingleCellExperiment(x), ...)
    }
)



# assays =======================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @export
setMethod(
    f = "assays",
    signature = signature("seurat"),
    definition = function(x, ...) {
        assays(as.SingleCellExperiment(x), ...)
    }
)



# colData ======================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
    f = "colData",
    signature = signature("seurat"),
    definition = function(x, ...) {
        colData(as.SingleCellExperiment(x), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(
    f = "colData<-",
    signature = signature(
        x = "seurat",
        value = "DataFrame"
    ),
    definition = function(x, value) {
        slot(x, "meta.data") <- as.data.frame(value)
        validObject(x)
        x
    }
)



# colnames =====================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics colnames
#' @export
setMethod(
    f = "colnames",
    signature = signature("seurat"),
    definition = function(x) {
        colnames(as.SingleCellExperiment(x))
    }
)



# counts =======================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics counts
#' @export
setMethod(
    f = "counts",
    signature = signature("seurat"),
    definition = function(object, ...) {
        counts(as.SingleCellExperiment(object), ...)
    }
)



# Gene2Symbol ==================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("seurat"),
    definition = function(object, ...) {
        Gene2Symbol(as(object, "SummarizedExperiment"), ...)
    }
)



# interestingGroups ============================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups
#' @export
setMethod(
    f = "interestingGroups",
    signature = signature("seurat"),
    definition = function(object, ...) {
        interestingGroups(as(object, "SummarizedExperiment"), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups<-
#' @export
setMethod(
    f = "interestingGroups<-",
    signature = signature(
        object = "seurat",
        value = "character"
    ),
    getMethod(
        f = "interestingGroups<-",
        signature = signature(
            object = "SummarizedExperiment",
            value = "character"
        ),
        where = asNamespace("basejump")
    )
)



# mapGenes =====================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump mapGenesToIDs
#' @export
setMethod(
    f = "mapGenesToIDs",
    signature = signature("seurat"),
    definition = function(object, genes, strict = TRUE) {
        mapGenesToIDs(
            object = as(object, "RangedSummarizedExperiment"),
            genes = genes,
            strict = strict
        )
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump mapGenesToRownames
#' @export
setMethod(
    f = "mapGenesToRownames",
    signature = signature("seurat"),
    definition = function(object, genes, strict = TRUE) {
        mapGenesToRownames(
            object = as(object, "RangedSummarizedExperiment"),
            genes = genes,
            strict = strict
        )
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump mapGenesToSymbols
#' @export
setMethod(
    f = "mapGenesToSymbols",
    signature = signature("seurat"),
    definition = function(object, genes, strict = TRUE) {
        mapGenesToSymbols(
            object = as(object, "RangedSummarizedExperiment"),
            genes = genes,
            strict = strict
        )
    }
)



# metadata =====================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    f = "metadata",
    signature = signature("seurat"),
    definition = function(x, ...) {
        stash <- .getSeuratStash(x, "metadata")
        if (!is.null(stash)) {
            return(stash)
        }
        metadata(as.SingleCellExperiment(x), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom S4Vectors metadata<-
#' @seealso `getMethod("metadata<-", "Annotated")`
#' @export
setMethod(
    f = "metadata<-",
    signature = signature(
        x = "seurat",
        value = "ANY"
    ),
    definition = function(x, value) {
        assert_is_list(value)
        if (!length(value)) {
            names(value) <- NULL
        }
        slot(x, "misc")[["metadata"]] <- value
        x
    }
)



# metrics ======================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("seurat"),
    definition = function(object, ...) {
        metrics(as(object, "SingleCellExperiment"), ...)
    }
)



# reducedDims ==================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDims
#' @export
setMethod(
    f = "reducedDims",
    signature = signature("seurat"),
    definition = function(x, ...) {
        reducedDims(as.SingleCellExperiment(x), ...)
    }
)



# rowData ======================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowData
#' @export
setMethod(
    f = "rowData",
    signature = signature("seurat"),
    definition = function(x, ...) {
        rowData(as(x, "SingleCellExperiment"), ...)
    }
)



# rownames =====================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics rownames
#' @export
setMethod(
    f = "rownames",
    signature = signature("seurat"),
    definition = function(x) {
        rownames(as.SingleCellExperiment(x))
    }
)



# rowRanges ====================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowRanges
#' @export
setMethod(
    f = "rowRanges",
    signature = signature("seurat"),
    definition = function(x) {
        sce <- as.SingleCellExperiment(x)
        # Default coercion method will return a GRangesList.
        gr <- rowRanges(sce)
        # Attempt to use stashed rowRanges, if defined.
        stash <- .getSeuratStash(x, "rowRanges")
        if (is(stash, "GRanges")) {
            assert_are_identical(names(gr), names(stash))
            mcols1 <- mcols(stash)
            mcols2 <- mcols(gr)
            mcols2 <- mcols2[
                ,
                setdiff(colnames(mcols2), colnames(mcols1)),
                drop = FALSE
            ]
            mcols <- cbind(mcols1, mcols2)
            gr <- stash
            mcols(gr) <- mcols
        }
        gr
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowRanges<-
#' @export
setMethod(
    f = "rowRanges<-",
    signature = signature("seurat"),
    definition = function(x, value) {
        assert_are_identical(rownames(x), names(value))
        misc <- slot(x, "misc")
        misc[["rowRanges"]] <- value
        slot(x, "misc") <- misc
        x
    }
)



# sampleData ===================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump sampleData
#' @export
setMethod(
    f = "sampleData",
    signature = signature("seurat"),
    getMethod(
        f = "sampleData",
        signature = signature("SingleCellExperiment"),
        where = asNamespace("basejump")
    )
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump sampleData<-
#' @export
setMethod(
    f = "sampleData<-",
    signature = signature(
        object = "seurat",
        value = "DataFrame"
    ),
    getMethod(
        f = "sampleData<-",
        signature = signature(
            object = "SingleCellExperiment",
            value = "DataFrame"
        ),
        where = asNamespace("basejump")
    )
)



# sampleNames ==================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump sampleNames
#' @export
setMethod(
    f = "sampleNames",
    signature = signature("seurat"),
    definition = function(object) {
        sampleNames(as(object, "SingleCellExperiment"))
    }
)



# weights ======================================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SingleCellExperiment weights
#' @export
setMethod(
    f = "weights",
    signature = signature("seurat"),
    definition = function(object) {
        slot(object, name = "misc")[["assays"]][["weights"]]
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SingleCellExperiment weights<-
#' @export
setMethod(
    f = "weights<-",
    signature = signature("seurat"),
    definition = function(object, value) {
        slot(object, name = "misc")[["assays"]][["weights"]] <- value
        object
    }
)

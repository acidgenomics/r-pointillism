#' Extend S4 Methods for `Seurat` Class
#'
#' Provide limited `SingleCellExperiment`-like method support.
#'
#' @name Seurat-SingleCellExperiment
#' @keywords internal
#'
#' @inheritParams basejump::params
#'
#' @return Match `SummarizedExperiment` method return.
NULL



# Internal =====================================================================
.getSeuratStash <- function(object, name) {
    assert(
        is(object, "Seurat"),
        isString(name)
    )

    misc <- slot(object, name = "misc")

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
#' @rdname Seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(
    f = "assay",
    signature = signature("Seurat"),
    definition = function(x, ...) {
        assay(as.SingleCellExperiment(x), ...)
    }
)



# assayNames ===================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assayNames
#' @export
setMethod(
    f = "assayNames",
    signature = signature("Seurat"),
    definition = function(x, ...) {
        assayNames(as.SingleCellExperiment(x), ...)
    }
)



# assays =======================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @export
setMethod(
    f = "assays",
    signature = signature("Seurat"),
    definition = function(x, ...) {
        assays(as.SingleCellExperiment(x), ...)
    }
)



# colData ======================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
    f = "colData",
    signature = signature("Seurat"),
    definition = function(x, ...) {
        colData(as.SingleCellExperiment(x), ...)
    }
)



#' @rdname Seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(
    f = "colData<-",
    signature = signature(
        x = "Seurat",
        value = "DataFrame"
    ),
    definition = function(x, value) {
        slot(x, "meta.data", check = TRUE) <- as.data.frame(value)
        validObject(x)
        x
    }
)



# colnames =====================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom BiocGenerics colnames
#' @export
setMethod(
    f = "colnames",
    signature = signature("Seurat"),
    definition = function(x) {
        colnames(as.SingleCellExperiment(x))
    }
)



# counts =======================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom BiocGenerics counts
#' @export
setMethod(
    f = "counts",
    signature = signature("Seurat"),
    definition = function(object, ...) {
        counts(as.SingleCellExperiment(object), ...)
    }
)



# Gene2Symbol ==================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("Seurat"),
    definition = function(object, ...) {
        Gene2Symbol(as(object, "SummarizedExperiment"), ...)
    }
)



# interestingGroups ============================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups
#' @export
setMethod(
    f = "interestingGroups",
    signature = signature("Seurat"),
    definition = function(object, ...) {
        interestingGroups(as(object, "SummarizedExperiment"), ...)
    }
)



#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups<-
#' @export
setMethod(
    f = "interestingGroups<-",
    signature = signature(
        object = "Seurat",
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
#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump mapGenesToIDs
#' @export
setMethod(
    f = "mapGenesToIDs",
    signature = signature("Seurat"),
    definition = function(object, genes, strict = TRUE) {
        mapGenesToIDs(
            object = as(object, "RangedSummarizedExperiment"),
            genes = genes,
            strict = strict
        )
    }
)



#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump mapGenesToRownames
#' @export
setMethod(
    f = "mapGenesToRownames",
    signature = signature("Seurat"),
    definition = function(object, genes, strict = TRUE) {
        mapGenesToRownames(
            object = as(object, "RangedSummarizedExperiment"),
            genes = genes,
            strict = strict
        )
    }
)



#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump mapGenesToSymbols
#' @export
setMethod(
    f = "mapGenesToSymbols",
    signature = signature("Seurat"),
    definition = function(object, genes, strict = TRUE) {
        mapGenesToSymbols(
            object = as(object, "RangedSummarizedExperiment"),
            genes = genes,
            strict = strict
        )
    }
)



# metadata =====================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    f = "metadata",
    signature = signature("Seurat"),
    definition = function(x, ...) {
        stash <- .getSeuratStash(x, "metadata")
        if (!is.null(stash)) {
            return(stash)
        }
        metadata(as.SingleCellExperiment(x), ...)
    }
)



#' @rdname Seurat-SingleCellExperiment
#' @importFrom S4Vectors metadata<-
#' @seealso `getMethod("metadata<-", "Annotated")`
#' @export
setMethod(
    f = "metadata<-",
    signature = signature(
        x = "Seurat",
        value = "ANY"
    ),
    definition = function(x, value) {
        assert(is.list(value))
        if (!length(value)) {
            names(value) <- NULL
        }
        x@misc[["metadata"]] <- value
        x
    }
)



# metrics ======================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("Seurat"),
    definition = function(object, ...) {
        metrics(as(object, "SingleCellExperiment"), ...)
    }
)



# reducedDims ==================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDims
#' @export
setMethod(
    f = "reducedDims",
    signature = signature("Seurat"),
    definition = function(x, ...) {
        reducedDims(as.SingleCellExperiment(x), ...)
    }
)



# rowData ======================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowData
#' @export
setMethod(
    f = "rowData",
    signature = signature("Seurat"),
    definition = function(x, ...) {
        rowData(as(x, "SingleCellExperiment"), ...)
    }
)



# rownames =====================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom BiocGenerics rownames
#' @export
setMethod(
    f = "rownames",
    signature = signature("Seurat"),
    definition = function(x) {
        rownames(as.SingleCellExperiment(x))
    }
)



# rowRanges ====================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowRanges
#' @export
setMethod(
    f = "rowRanges",
    signature = signature("Seurat"),
    definition = function(x) {
        sce <- as.SingleCellExperiment(x)
        # Default coercion method will return a GRangesList.
        gr <- rowRanges(sce)
        # Attempt to use stashed rowRanges, if defined.
        stash <- .getSeuratStash(x, "rowRanges")
        if (is(stash, "GRanges")) {
            assert(identical(length(gr), length(stash)))
            # Handle situation where we've changed from gene IDs to gene names.
            if (!identical(names(gr), names(stash))) {
                names(stash) <- names(gr)
            }
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



#' @rdname Seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowRanges<-
#' @export
setMethod(
    f = "rowRanges<-",
    signature = signature("Seurat"),
    definition = function(x, value) {
        assert(identical(rownames(x), names(value)))
        x@misc[["rowRanges"]] <- value
        x
    }
)



# sampleData ===================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump sampleData
#' @export
setMethod(
    f = "sampleData",
    signature = signature("Seurat"),
    getMethod(
        f = "sampleData",
        signature = signature("SingleCellExperiment"),
        where = asNamespace("basejump")
    )
)



#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump sampleData<-
#' @export
setMethod(
    f = "sampleData<-",
    signature = signature(
        object = "Seurat",
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
#' @rdname Seurat-SingleCellExperiment
#' @importFrom basejump sampleNames
#' @export
setMethod(
    f = "sampleNames",
    signature = signature("Seurat"),
    definition = function(object) {
        sampleNames(as(object, "SingleCellExperiment"))
    }
)



# weights ======================================================================
#' @rdname Seurat-SingleCellExperiment
#' @importFrom SingleCellExperiment weights
#' @export
setMethod(
    f = "weights",
    signature = signature("Seurat"),
    definition = function(object) {
        slot(object, name = "misc")[["assays"]][["weights"]]
    }
)



#' @rdname Seurat-SingleCellExperiment
#' @importFrom SingleCellExperiment weights<-
#' @export
setMethod(
    f = "weights<-",
    signature = signature("Seurat"),
    definition = function(object, value) {
        object@misc[["assays"]][["weights"]] <- value
        object
    }
)

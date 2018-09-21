# FIXME
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function 'mapGenesToSymbols' for signature '"seurat"'



#' Extend S4 Methods for `seurat` Class
#'
#' Provide limited `SingleCellExperiment`-like method support.
#'
#' @name seurat-SingleCellExperiment
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return Match `SummarizedExperiment` method return.
NULL



.getSeuratStash <- function(object, name) {
    stopifnot(is(object, "seurat"))
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



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(
    f = "assay",
    signature = signature("seurat"),
    definition = function(x, ...) {
        assay(.as.SingleCellExperiment.seurat(x), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assayNames
#' @export
setMethod(
    f = "assayNames",
    signature = signature("seurat"),
    definition = function(x, ...) {
        assayNames(.as.SingleCellExperiment.seurat(x), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @export
setMethod(
    f = "assays",
    signature = signature("seurat"),
    definition = function(x, ...) {
        assays(.as.SingleCellExperiment.seurat(x), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
    f = "colData",
    signature = signature("seurat"),
    definition = function(x, ...) {
        colData(.as.SingleCellExperiment.seurat(x), ...)
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



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics colnames
#' @export
setMethod(
    f = "colnames",
    signature = signature("seurat"),
    definition = function(x) {
        colnames(.as.SingleCellExperiment.seurat(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics counts
#' @export
setMethod(
    f = "counts",
    signature = signature("seurat"),
    definition = function(object, ...) {
        counts(.as.SingleCellExperiment.seurat(object), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump gene2symbol
#' @export
setMethod(
    f = "gene2symbol",
    signature = signature("seurat"),
    definition = function(object) {
        gene2symbol(as(object, "SummarizedExperiment"))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups
#' @export
setMethod(
    f = "interestingGroups",
    signature = signature("seurat"),
    definition = function(object) {
        interestingGroups(as(object, "SummarizedExperiment"))
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
        metadata(.as.SingleCellExperiment.seurat(x), ...)
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



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("seurat"),
    definition = function(object, ...) {
        metrics(as(object, "SingleCellExperiment"))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDims
#' @export
setMethod(
    f = "reducedDims",
    signature = signature("seurat"),
    definition = function(x) {
        reducedDims(.as.SingleCellExperiment.seurat(x))
    }
)



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



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics rownames
#' @export
setMethod(
    f = "rownames",
    signature = signature("seurat"),
    definition = function(x) {
        rownames(.as.SingleCellExperiment.seurat(x))
    }
)



# Note that there may be a names issue when handling GRanges with NA values
# in either the geneID or geneName columns
#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowRanges
#' @export
setMethod(
    f = "rowRanges",
    signature = signature("seurat"),
    definition = function(x) {
        gr <- rowRanges(.as.SingleCellExperiment.seurat(x))
        # Attempt to use stashed rowRanges, if defined.
        stash <- .getSeuratStash(x, "rowRanges")
        if (is(stash, "GRanges")) {
            assert_is_subset(c("geneID", "geneName"), colnames(mcols(stash)))
            # Check to see if we're using IDs or symbols
            if (any(names(gr) %in% mcols(stash)[["geneID"]])) {
                namesCol <- "geneID"
            } else if (any(names(gr) %in% mcols(stash)[["geneName"]])) {
                namesCol <- "geneName"
            } else {
                stop("Failed to match identifiers to rownames")
            }
            names(stash) <- make.unique(as.character(mcols(stash)[[namesCol]]))
            assert_are_identical(names(gr), names(stash))
            assert_are_disjoint_sets(
                x = colnames(mcols(gr)),
                y = colnames(mcols(stash))
            )
            mcols(stash) <- cbind(mcols(stash), mcols(gr))
            gr <- stash
        }
        gr
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioSingleCell sampleData
#' @export
setMethod(
    f = "sampleData",
    signature = signature("seurat"),
    getMethod(
        f = "sampleData",
        signature = signature("SingleCellExperiment"),
        where = asNamespace("bcbioSingleCell")
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
        where = asNamespace("bcbioSingleCell")
    )
)



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

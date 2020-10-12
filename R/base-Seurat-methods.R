## Consider adding `rowData<-` assignment method.



#' Extend base S4 methods for `Seurat` class
#'
#' Provide limited `SingleCellExperiment`-like method support.
#'
#' @name base-Seurat
#' @keywords internal
#' @note Updated 2019-08-05.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @seealso
#' - `Seurat::GetAssayData()`.
#'
#' @return Varies, depending on the generic.
NULL



## Updated 2019-08-02.
`Gene2Symbol,Seurat` <-  # nolint
    function(object, ...) {
        Gene2Symbol(
            object = as(object, "SummarizedExperiment"),
            ...
        )
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("Seurat"),
    definition = `Gene2Symbol,Seurat`
)



## Updated 2019-08-05.
`assay,Seurat` <-  # nolint
    function(x, ...) {
        x <- as.SingleCellExperiment(x)
        assay(x, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "assay",
    signature = signature("Seurat"),
    definition = `assay,Seurat`
)



## Updated 2019-08-05.
`assayNames,Seurat` <-  # nolint
    function(x, ...) {
        x <- as.SingleCellExperiment(x)
        assayNames(x, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "assayNames",
    signature = signature("Seurat"),
    definition = `assayNames,Seurat`
)



## Updated 2019-08-05.
`assays,Seurat` <-  # nolint
    function(x, ...) {
        x <- as.SingleCellExperiment(x)
        assays(x, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "assays",
    signature = signature("Seurat"),
    definition = `assays,Seurat`
)



## Updated 2019-08-05.
`colData,Seurat` <-  # nolint
    function(x, ...) {
        x <- as.SingleCellExperiment(x)
        colData(x, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "colData",
    signature = signature("Seurat"),
    definition = `colData,Seurat`
)



## Updated 2019-08-05.
`colData<-,Seurat,DataFrame` <-  # nolint
    function(x, value) {
        value <- as.data.frame(value)
        slot(x, "meta.data", check = TRUE) <- value
        validObject(x)
        x
    }



#' @rdname base-Seurat
#' @export
setReplaceMethod(
    f = "colData",
    signature = signature(
        x = "Seurat",
        value = "DataFrame"
    ),
    definition = `colData<-,Seurat,DataFrame`
)



## Upated 2019-08-06.
`counts,Seurat` <-  # nolint
    function(object, assay = NULL) {
        object <- as.SingleCellExperiment(object, assay = assay)
        counts(object)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "counts",
    signature = signature("Seurat"),
    definition = `counts,Seurat`
)



## Updated 2019-08-05.
`interestingGroups,Seurat` <-  # nolint
    function(object, ...) {
        object <- as(object, "SingleCellExperiment")
        interestingGroups(object, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "interestingGroups",
    signature = signature("Seurat"),
    definition = `interestingGroups,Seurat`
)



## Updated 2019-08-05.
`interestingGroups<-,Seurat,character` <-  # nolint
    getMethod(
        f = "interestingGroups<-",
        signature = signature(
            object = "SummarizedExperiment",
            value = "character"
        ),
        where = asNamespace("basejump")
    )



#' @rdname base-Seurat
#' @export
setReplaceMethod(
    f = "interestingGroups",
    signature = signature(
        object = "Seurat",
        value = "character"
    ),
    definition = `interestingGroups<-,Seurat,character`
)



## Updated 2020-01-30.
`logcounts,Seurat` <-  # nolint
    function(object, assay = NULL) {
        norm <- .seuratNormalizationMethod(object, assay = assay)
        if (norm != "LogNormalize") {
            object <- NormalizeData(
                object = object,
                normalization.method = "LogNormalize",
                verbose = TRUE
            )
        }
        GetAssayData(object = object, assay = assay, slot = "data")
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "logcounts",
    signature = signature("Seurat"),
    definition = `logcounts,Seurat`
)



## Updated 2019-08-05.
`mapGenesToIDs,Seurat` <-  # nolint
    function(object, genes, strict = TRUE) {
        object <- as(object, "SingleCellExperiment")
        mapGenesToIDs(
            object = object,
            genes = genes,
            strict = strict
        )
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "mapGenesToIDs",
    signature = signature("Seurat"),
    definition = `mapGenesToIDs,Seurat`
)



## Updated 2019-08-05.
`mapGenesToRownames,Seurat` <-  # nolint
    function(object, genes, strict = TRUE) {
        object <- as(object, "SingleCellExperiment")
        mapGenesToRownames(
            object = object,
            genes = genes,
            strict = strict
        )
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "mapGenesToRownames",
    signature = signature("Seurat"),
    definition = `mapGenesToRownames,Seurat`
)



## Updated 2019-08-05.
`mapGenesToSymbols,Seurat` <-  # nolint
    function(object, genes, strict = TRUE) {
        object <- as(object, "SingleCellExperiment")
        mapGenesToSymbols(
            object = object,
            genes = genes,
            strict = strict
        )
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "mapGenesToSymbols",
    signature = signature("Seurat"),
    definition = `mapGenesToSymbols,Seurat`
)



## Updated 2019-08-05.
`metadata,Seurat` <-  # nolint
    function(x, ...) {
        stash <- .getSeuratStash(x, "metadata")
        if (!is.null(stash)) {
            return(stash)
        }
        x <- as.SingleCellExperiment(x)
        metadata(x, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "metadata",
    signature = signature("Seurat"),
    definition = `metadata,Seurat`
)



## Updated 2019-08-05.
`metadata<-,Seurat,ANY` <-  # nolint
    function(x, value) {
        assert(is.list(value))
        if (!length(value)) {
            names(value) <- NULL
        }
        x@misc[["metadata"]] <- value
        x
    }



#' @rdname base-Seurat
#' @export
setReplaceMethod(
    f = "metadata",
    signature = signature(
        x = "Seurat",
        value = "ANY"
    ),
    definition = `metadata<-,Seurat,ANY`
)



## Updated 2019-08-05.
`metrics,Seurat` <-  # nolint
    function(object, ...) {
        object <- as(object, "SingleCellExperiment")
        metrics(object, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "metrics",
    signature = signature("Seurat"),
    definition = `metrics,Seurat`
)



## Updated 2020-01-30.
`normcounts,Seurat` <-  # nolint
    function(object, assay = NULL) {
        ## Check for pre-calculated relative counts (not typical).
        method <- .seuratNormalizationMethod(object, assay = assay)
        scaleFactor <- .seuratScaleFactor(object, assay = assay)
        if (!(method == "RC" && scaleFactor == 1L)) {
            object <- NormalizeData(
                object = object,
                normalization.method = "RC",
                scale.factor = 1L,
                verbose = TRUE
            )
        }
        GetAssayData(object = object, assay = assay, slot = "data")
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "normcounts",
    signature = signature("Seurat"),
    definition = `normcounts,Seurat`
)



## Updated 2019-08-06.
`organism,Seurat` <-  # nolint
    function(object) {
        object <- as(object, "SingleCellExperiment")
        organism(object)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "organism",
    signature = signature("Seurat"),
    definition = `organism,Seurat`
)



## Updated 2019-08-06.
`organism<-,Seurat,character` <-  # nolint
    function(object, value) {
        metadata(object)[["organism"]] <- value
        object
    }



#' @rdname base-Seurat
#' @export
setReplaceMethod(
    f = "organism",
    signature = "Seurat",
    definition = `organism<-,Seurat,character`
)



## Updated 2019-08-06.
`reducedDim,Seurat` <-  # nolint
    function(x, type = 1L, withDimnames = TRUE) {
        x <- as.SingleCellExperiment(x)
        reducedDim(x = x, type = type, withDimnames = withDimnames)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "reducedDim",
    signature = signature("Seurat"),
    definition = `reducedDim,Seurat`
)



## Updated 2019-08-05.
`reducedDimNames,Seurat` <-  # nolint
    function(x) {
        x <- as.SingleCellExperiment(x)
        reducedDimNames(x)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "reducedDimNames",
    signature = signature("Seurat"),
    definition = `reducedDimNames,Seurat`
)



## Updated 2019-08-05.
`reducedDims,Seurat` <-  # nolint
    function(x, ...) {
        x <- as.SingleCellExperiment(x)
        reducedDims(x, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "reducedDims",
    signature = signature("Seurat"),
    definition = `reducedDims,Seurat`
)



## Updated 2019-08-05.
`rowData,Seurat` <-  # nolint
    function(x, ...) {
        x <- as(x, "SingleCellExperiment")
        rowData(x, ...)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "rowData",
    signature = signature("Seurat"),
    definition = `rowData,Seurat`
)



## Updated 2020-02-20.
`rowRanges,Seurat` <-  # nolint
    function(x) {
        sce <- as.SingleCellExperiment(x)
        gr <- rowRanges(sce)
        ## Attempt to use stashed rowRanges, if defined.
        stash <- .getSeuratStash(x, "rowRanges")
        if (!is(stash, "GRanges")) {
            return(gr)
        }
        ## Handle situation where we've changed from gene IDs to gene names
        ## using `convertGenesToSymbols`.
        if (
            identical(length(gr), length(stash)) &&
            !identical(names(gr), names(stash))
        ) {
            names(stash) <- names(gr)
        }
        assert(areIntersectingSets(names(stash), names(gr)))
        stash <- stash[names(gr)]
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
        gr
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "rowRanges",
    signature = signature("Seurat"),
    definition = `rowRanges,Seurat`
)



## Updated 2019-08-05.
`rowRanges<-,Seurat,ANY` <-  # nolint
    function(x, value) {
        assert(identical(rownames(x), names(value)))
        x@misc[["rowRanges"]] <- value
        x
    }



#' @rdname base-Seurat
#' @export
setReplaceMethod(
    f = "rowRanges",
    signature = signature(
        x = "Seurat",
        value = "ANY"
    ),
    definition = `rowRanges<-,Seurat,ANY`
)



## Updated 2019-08-05.
`sampleData,Seurat` <-  # nolint
    getMethod(
        f = "sampleData",
        signature = signature("SingleCellExperiment"),
        where = asNamespace("basejump")
    )



#' @rdname base-Seurat
#' @export
setMethod(
    f = "sampleData",
    signature = signature("Seurat"),
    definition = `sampleData,Seurat`
)



## Updated 2019-08-05.
`sampleData<-,Seurat,DataFrame` <-  # nolint
    getMethod(
        f = "sampleData<-",
        signature = signature(
            object = "SingleCellExperiment",
            value = "DataFrame"
        ),
        where = asNamespace("basejump")
    )



#' @rdname base-Seurat
#' @export
setReplaceMethod(
    f = "sampleData",
    signature = signature(
        object = "Seurat",
        value = "DataFrame"
    ),
    definition = `sampleData<-,Seurat,DataFrame`
)



## Updated 2019-08-05.
`sampleNames,Seurat` <-  # nolint
    function(object) {
        object <- as(object, "SingleCellExperiment")
        sampleNames(object)
    }



#' @rdname base-Seurat
#' @export
setMethod(
    f = "sampleNames",
    signature = signature("Seurat"),
    definition = `sampleNames,Seurat`
)

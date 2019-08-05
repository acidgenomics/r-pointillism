## Consider adding `rowData<-` assignment support.



#' Extend base S4 methods for `Seurat` class
#'
#' Provide limited `SingleCellExperiment`-like method support.
#'
#' @name Seurat-SingleCellExperiment
#' @keywords internal
#' @note Updated 2019-08-05.
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SingleCellExperiment reducedDims reducedDimNames
#' @importFrom SummarizedExperiment assay assayNames assays colData colData<-
#'   rowData rowRanges rowRanges<-
#' @importFrom basejump interestingGroups interestingGroups<- mapGenesToIDs
#'   mapGenesToRownames mapGenesToSymbols metrics reducedDims sampleData
#'   sampleData<- sampleNames
#'
#' @inheritParams basejump::params
#'
#' @seealso
#' - `Seurat::GetAssayData()`.
#'
#' @return Match `SingleCellExperiment` method return.
NULL



## Updated 2019-08-05.
`assay,Seurat` <-  # nolint
    function(x, ...) {
        x <- as.SingleCellExperiment(x)
        assay(x, ...)
    }



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
#' @export
setReplaceMethod(
    f = "colData",
    signature = signature(
        x = "Seurat",
        value = "DataFrame"
    ),
    definition = `colData<-,Seurat,DataFrame`
)



## Updated 2019-08-05.
`interestingGroups,Seurat` <-  # nolint
    function(object, ...) {
        object <- as(object, "SingleCellExperiment")
        interestingGroups(object, ...)
    }



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
#' @export
setReplaceMethod(
    f = "interestingGroups",
    signature = signature(
        object = "Seurat",
        value = "character"
    ),
    definition = `interestingGroups<-,Seurat,character`
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
#' @export
setMethod(
    f = "metrics",
    signature = signature("Seurat"),
    definition = `metrics,Seurat`
)



## Updated 2019-08-05.
`reducedDimNames,Seurat` <-  # nolint
    function(x) {
        x <- as.SingleCellExperiment(x)
        reducedDimNames(x)
    }



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
#' @export
setMethod(
    f = "rowData",
    signature = signature("Seurat"),
    definition = `rowData,Seurat`
)



## Updated 2019-08-05.
`rowRanges,Seurat` <-  # nolint
    function(x) {
        sce <- as.SingleCellExperiment(x)
        ## Default coercion method will return a GRangesList.
        gr <- rowRanges(sce)
        ## Attempt to use stashed rowRanges, if defined.
        stash <- .getSeuratStash(x, "rowRanges")
        if (is(stash, "GRanges")) {
            assert(identical(length(gr), length(stash)))
            ## Handle situation where we've changed from gene IDs to gene names.
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
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



#' @rdname Seurat-SingleCellExperiment
#' @export
setMethod(
    f = "sampleNames",
    signature = signature("Seurat"),
    definition = `sampleNames,Seurat`
)

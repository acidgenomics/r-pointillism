#' Extend base S4 methods for `Seurat` class
#'
#' Provide limited `SingleCellExperiment`-like method support.
#'
#' @name base-Seurat
#' @keywords internal
#' @note Updated 2023-08-16.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @seealso
#' - `SeuratObject::LayerData()`.
#'
#' @return Varies, depending on the generic.
NULL



## Updated 2019-08-02.
`GeneToSymbol,Seurat` <- # nolint
    function(object, ...) {
        GeneToSymbol(
            object = as(object, "SummarizedExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "GeneToSymbol",
    signature = signature(object = "Seurat"),
    definition = `GeneToSymbol,Seurat`
)



## Updated 2021-03-03.
`assay,Seurat` <- # nolint
    function(x, ...) {
        assay(
            x = as(x, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "assay",
    signature = signature(x = "Seurat"),
    definition = `assay,Seurat`
)



## Updated 2021-03-03.
`assayNames,Seurat` <- # nolint
    function(x, ...) {
        assayNames(
            x = as(x, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "assayNames",
    signature = signature(x = "Seurat"),
    definition = `assayNames,Seurat`
)



## Updated 2021-03-03.
`assays,Seurat` <- # nolint
    function(x, ...) {
        assays(
            x = as(x, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "assays",
    signature = signature(x = "Seurat"),
    definition = `assays,Seurat`
)



## Updated 2021-09-13.
`cellCountsPerCluster,Seurat` <- # nolint
    function(object, ...) {
        cellCountsPerCluster(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "cellCountsPerCluster",
    signature = signature(object = "Seurat"),
    definition = `cellCountsPerCluster,Seurat`
)



## Updated 2019-08-02.
`clusters,Seurat` <- # nolint
    function(object) {
        validObject(object)
        Idents(object)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "clusters",
    signature = signature(object = "Seurat"),
    definition = `clusters,Seurat`
)



## Updated 2023-08-16.
`colData,Seurat` <- # nolint
    function(x) {
        assert(validObject(x))
        df <- slot(object = x, name = "meta.data")
        assert(
            is.data.frame(df),
            identical(colnames(x), rownames(df))
        )
        df <- as(df, "DFrame")
        df
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "colData",
    signature = signature(x = "Seurat"),
    definition = `colData,Seurat`
)



## Updated 2023-08-16.
`colData<-,Seurat,DFrame` <- # nolint
    function(x, value) {
        value <- as.data.frame(value)
        slot(x, name = "meta.data", check = TRUE) <- value # nolint
        assert(validObject(x))
        x
    }

#' @rdname base-Seurat
#' @export
setReplaceMethod(
    f = "colData",
    signature = signature(
        x = "Seurat",
        value = "DFrame"
    ),
    definition = `colData<-,Seurat,DFrame`
)



## Upated 2023-08-16.
`counts,Seurat` <- # nolint
    function(object, assay = NULL) {
        LayerData(object, layer = "counts", assay = assay)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "counts",
    signature = signature(object = "Seurat"),
    definition = `counts,Seurat`
)



## Updated 2023-08-16.
`cpm,Seurat` <- # nolint
    function(object, assay = NULL) {
        ## Check for pre-calculated CPM (not typical).
        method <- .seuratNormalizationMethod(object, assay = assay)
        scaleFactor <- .seuratScaleFactor(object, assay = assay)
        if (!(method == "RC" && scaleFactor == 1e6L)) {
            alert(
                "Generating CPM with {.pkg Seurat}::{.fun NormalizeData}."
            )
            object <- NormalizeData(
                object = object,
                assay = assay,
                normalization.method = "RC",
                scale.factor = 1e6L,
                verbose = TRUE
            )
        }
        LayerData(object = object, layer = "data", assay = assay)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "cpm",
    signature = signature(object = "Seurat"),
    definition = `cpm,Seurat`
)



## Updated 2020-02-21.
`diffExp,Seurat` <- # nolint
    function(object, ...) {
        diffExp(object = as(object, "SingleCellExperiment"), ...)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "diffExp",
    signature = signature(object = "Seurat"),
    definition = `diffExp,Seurat`
)



## Updated 2021-09-13.
`diffExpPerCluster,Seurat` <- # nolint
    function(object, ...) {
        diffExpPerCluster(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "diffExpPerCluster",
    signature = signature(object = "Seurat"),
    definition = `diffExpPerCluster,Seurat`
)



## Updated 2021-09-13.
`findMarkers,Seurat` <- # nolint
    function(object, ...) {
        findMarkers(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "findMarkers",
    signature = signature(object = "Seurat"),
    definition = `findMarkers,Seurat`
)



## Updated 2021-03-03.
`interestingGroups,Seurat` <- # nolint
    function(object, ...) {
        interestingGroups(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "interestingGroups",
    signature = signature(object = "Seurat"),
    definition = `interestingGroups,Seurat`
)



## Updated 2021-03-02.
`interestingGroups<-,Seurat,character` <- # nolint
    getMethod(
        f = "interestingGroups<-",
        signature = signature(
            object = "SummarizedExperiment",
            value = "character"
        ),
        where = asNamespace("AcidExperiment")
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



## Updated 2023-08-16.
`logcounts,Seurat` <- # nolint
    function(object, assay = NULL) {
        norm <- .seuratNormalizationMethod(object, assay = assay)
        if (norm != "LogNormalize") {
            object <- NormalizeData(
                object = object,
                normalization.method = "LogNormalize",
                verbose = TRUE
            )
        }
        LayerData(object = object, layer = "data", assay = assay)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "logcounts",
    signature = signature(object = "Seurat"),
    definition = `logcounts,Seurat`
)



## Updated 2021-03-03.
`mapGenesToIDs,Seurat` <- # nolint
    function(object, genes, strict = TRUE) {
        mapGenesToIDs(
            object = as(object, "SingleCellExperiment"),
            genes = genes,
            strict = strict
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "mapGenesToIDs",
    signature = signature(object = "Seurat"),
    definition = `mapGenesToIDs,Seurat`
)



## Updated 2021-03-03.
`mapGenesToRownames,Seurat` <- # nolint
    function(object, genes, strict = TRUE) {
        mapGenesToRownames(
            object = as(object, "SingleCellExperiment"),
            genes = genes,
            strict = strict
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "mapGenesToRownames",
    signature = signature(object = "Seurat"),
    definition = `mapGenesToRownames,Seurat`
)



## Updated 2021-03-03.
`mapGenesToSymbols,Seurat` <- # nolint
    function(object, genes, strict = TRUE) {
        mapGenesToSymbols(
            object = as(object, "SingleCellExperiment"),
            genes = genes,
            strict = strict
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "mapGenesToSymbols",
    signature = signature(object = "Seurat"),
    definition = `mapGenesToSymbols,Seurat`
)



## Updated 2022-05-11.
`metadata,Seurat` <- # nolint
    function(x, ...) {
        stash <- .getSeuratStash(x, "metadata")
        if (!is.null(stash)) {
            return(stash)
        }
        x <- Seurat::as.SingleCellExperiment(x)
        metadata(x, ...)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "metadata",
    signature = signature(x = "Seurat"),
    definition = `metadata,Seurat`
)



## Updated 2019-08-05.
`metadata<-,Seurat,ANY` <- # nolint
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



## Updated 2021-03-03.
`metrics,Seurat` <- # nolint
    function(object, ...) {
        metrics(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "metrics",
    signature = signature(object = "Seurat"),
    definition = `metrics,Seurat`
)



## Updated 2020-01-30.
`normalize,Seurat` <- # nolint
    function(object) {
        alert(sprintf(
            "Normalizing with {.pkg %s}::{.fun %s}.",
            "Seurat", "NormalizeData"
        ))
        NormalizeData(object = object, verbose = TRUE)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "normalize",
    signature = signature(object = "Seurat"),
    definition = `normalize,Seurat`
)



## Updated 2023-08-16.
`normcounts,Seurat` <- # nolint
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
        LayerData(object = object, layer = "data", assay = assay)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "normcounts",
    signature = signature(object = "Seurat"),
    definition = `normcounts,Seurat`
)



## Updated 2021-03-03.
`organism,Seurat` <- # nolint
    function(object) {
        organism(object = as(object, "SingleCellExperiment"))
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "organism",
    signature = signature(object = "Seurat"),
    definition = `organism,Seurat`
)



## Updated 2019-08-06.
`organism<-,Seurat,character` <- # nolint
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



## Updated 2021-09-13.
`plotCellCountsPerCluster,Seurat` <- # nolint
    function(object, ...) {
        plotCellCountsPerCluster(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotCellCountsPerCluster",
    signature = signature(object = "Seurat"),
    definition = `plotCellCountsPerCluster,Seurat`
)



## Updated 2021-09-13.
`plotCellTypesPerCluster,Seurat,KnownMarkers` <- # nolint
    function(object, markers, ...) {
        plotCellTypesPerCluster(
            object = as(object, "SingleCellExperiment"),
            markers = markers,
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotCellTypesPerCluster",
    signature = signature(
        object = "Seurat",
        markers = "KnownMarkers"
    ),
    definition = `plotCellTypesPerCluster,Seurat,KnownMarkers`
)



## Updated 2021-09-13.
`plotCounts,Seurat` <- # nolint
    function(object, ...) {
        plotCounts(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotCounts",
    signature = signature(object = "Seurat"),
    definition = `plotCounts,Seurat`
)



## Updated 2021-09-13.
`plotDots,Seurat` <- # nolint
    function(object, ...) {
        plotDots(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotDots",
    signature = signature(object = "Seurat"),
    definition = `plotDots,Seurat`
)



## Updated 2021-09-13.
`plotFeature,Seurat` <- # nolint
    function(object, ...) {
        plotFeature(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotFeature",
    signature = signature(object = "Seurat"),
    definition = `plotFeature,Seurat`
)



## Updated 2021-09-13.
`plotKnownMarkers,Seurat,KnownMarkers` <- # nolint
    function(object, markers, ...) {
        plotKnownMarkers(
            object = as(object, "SingleCellExperiment"),
            markers = markers,
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotKnownMarkers",
    signature = signature(
        object = "Seurat",
        markers = "KnownMarkers"
    ),
    definition = `plotKnownMarkers,Seurat,KnownMarkers`
)



## Updated 2021-09-13.
`plotMarker,Seurat` <- # nolint
    function(object, ...) {
        plotMarker(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotMarker",
    signature = signature(object = "Seurat"),
    definition = `plotMarker,Seurat`
)



## Updated 2021-09-13.
`plotPCA,Seurat` <- # nolint
    function(object, ...) {
        plotPCA(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotPCA",
    signature = signature(object = "Seurat"),
    definition = `plotPCA,Seurat`
)



## Updated 2021-03-02.
`plotReducedDim,Seurat` <- # nolint
    function(object, ...) {
        validObject(object)
        idents <- .seuratWhichIdents(object)
        dl(c("idents" = idents))
        plotReducedDim(object = as(object, "SingleCellExperiment"), ...)
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotReducedDim",
    signature = signature(object = "Seurat"),
    definition = `plotReducedDim,Seurat`
)



## Updated 2021-09-13.
`plotStackedBarPlot,Seurat` <- # nolint
    function(object, ...) {
        plotStackedBarPlot(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotStackedBarPlot",
    signature = signature(object = "Seurat"),
    definition = `plotStackedBarPlot,Seurat`
)



## Updated 2021-09-13.
`plotTSNE,Seurat` <- # nolint
    function(object, ...) {
        plotTSNE(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotTSNE",
    signature = signature(object = "Seurat"),
    definition = `plotTSNE,Seurat`
)



## Updated 2021-09-13.
`plotViolin,Seurat` <- # nolint
    function(object, ...) {
        plotViolin(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotViolin",
    signature = signature(object = "Seurat"),
    definition = `plotViolin,Seurat`
)



## Updated 2021-09-13.
`plotUMAP,Seurat` <- # nolint
    function(object, ...) {
        plotUMAP(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "plotUMAP",
    signature = signature(object = "Seurat"),
    definition = `plotUMAP,Seurat`
)



## Updated 2021-03-03.
`reducedDim,Seurat` <- # nolint
    function(x, type = 1L, withDimnames = TRUE) {
        reducedDim(
            x = as(x, "SingleCellExperiment"),
            type = type,
            withDimnames = withDimnames
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "reducedDim",
    signature = signature(x = "Seurat"),
    definition = `reducedDim,Seurat`
)



## Updated 2021-03-03.
`reducedDimNames,Seurat` <- # nolint
    function(x) {
        reducedDimNames(x = as(x, "SingleCellExperiment"))
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "reducedDimNames",
    signature = signature(x = "Seurat"),
    definition = `reducedDimNames,Seurat`
)



## Updated 2021-03-03.
`reducedDims,Seurat` <- # nolint
    function(x, ...) {
        reducedDims(
            x = as(x, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "reducedDims",
    signature = signature(x = "Seurat"),
    definition = `reducedDims,Seurat`
)



## Updated 2021-03-03.
`rowData,Seurat` <- # nolint
    function(x, ...) {
        rowData(
            x = as(x, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "rowData",
    signature = signature(x = "Seurat"),
    definition = `rowData,Seurat`
)



## Updated 2022-05-11.
`rowRanges,Seurat` <- # nolint
    function(x) {
        sce <- Seurat::as.SingleCellExperiment(x)
        gr <- rowRanges(sce)
        assert(
            is(gr, "GenomicRangesList"),
            all(
                ## FIXME Can we replace with `lengths` here?
                vapply(
                    X = gr,
                    FUN = length,
                    FUN.VALUE = integer(1L)
                ) == 0L
            )
        )
        ## Attempt to use stashed rowRanges, if defined.
        ## Otherwise, return the empty placeholder.
        stash <- .getSeuratStash(x, "rowRanges")
        if (!is(stash, "GenomicRanges")) {
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
        colnames(mcols) <- camelCase(colnames(mcols), strict = TRUE)
        gr <- stash
        mcols(gr) <- mcols
        gr
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "rowRanges",
    signature = signature(x = "Seurat"),
    definition = `rowRanges,Seurat`
)



## Updated 2019-08-05.
`rowRanges<-,Seurat,ANY` <- # nolint
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



## Updated 2021-03-03.
`sampleData,Seurat` <- # nolint
    function(object, ...) {
        sampleData(
            object = as(object, "SingleCellExperiment"),
            ...
        )
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "sampleData",
    signature = signature(object = "Seurat"),
    definition = `sampleData,Seurat`
)



## Updated 2023-08-16.
`sampleData<-,Seurat,DFrame` <- # nolint
    getMethod(
        f = "sampleData<-",
        signature = signature(
            object = "SingleCellExperiment",
            value = "DFrame"
        ),
        where = asNamespace("AcidSingleCell")
    )

#' @rdname base-Seurat
#' @export
setReplaceMethod(
    f = "sampleData",
    signature = signature(
        object = "Seurat",
        value = "DFrame"
    ),
    definition = `sampleData<-,Seurat,DFrame`
)



## Updated 2021-03-03.
`sampleNames,Seurat` <- # nolint
    function(object) {
        sampleNames(object = as(object, "SingleCellExperiment"))
    }

#' @rdname base-Seurat
#' @export
setMethod(
    f = "sampleNames",
    signature = signature(object = "Seurat"),
    definition = `sampleNames,Seurat`
)

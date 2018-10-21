# TODO Improve the documentation here.



#' @inherit methods::as
#' @name as
#' @aliases coerce
#' @importFrom methods coerce
#' @exportMethod coerce
#'
#' @seealso
#' - [Seurat::CreateSeuratObject()].
#' - [Seurat::Convert()].
#'
#' @examples
#' data(seurat_small)
#' data(sce_small, package = "basejump")
#'
#' ## SingleCellExperiment to seurat ====
#' object <- sce_small
#' print(object)
#' x <- as(object, "seurat")
#' class(x)
#' print(x)
#'
#' ## seurat to SingleCellExperiment ====
#' x <- as(seurat_small, "SingleCellExperiment")
#' print(x)
NULL



# SingleCellExperiment to seurat ===============================================
as.seurat <- function(from) {
    UseMethod("as.seurat")
}



as.seurat.SingleCellExperiment <- function(from) {
    # Create the seurat object.
    args <- list(
        raw.data = counts(from),
        meta.data = as.data.frame(colData(from)),
        # Already applied filtering cutoffs for cells and genes.
        min.cells = 0L,
        min.genes = 0L,
        # Using the default for UMI datasets.
        is.expr = 0L,
        project = "pointillism"
    )
    to <- do.call(what = CreateSeuratObject, args = args)

    # Check that the dimensions match exactly.
    assert_are_identical(
        x = dim(from),
        y = dim(slot(to, "raw.data"))
    )

    # Keep extra assays, if defined (e.g. weights).
    assayNames <- setdiff(
        x = assayNames(from),
        y = c(
            "counts",
            "logcounts",
            "scaleData"  # from Seurat-to-SCE.
        )
    )
    assays <- assays(from)[assayNames]

    # rowRanges.
    rowRanges <- rowRanges(from)

    # metadata.
    metadata <- metadata(from)
    metadata[["sessionInfo"]] <- session_info()

    misc <- list(
        assays = assays,
        rowRanges = rowRanges,
        metadata = metadata
    )
    misc <- Filter(Negate(is.null), misc)
    slot(to, name = "misc") <- misc

    to
}



#' @rdname as
#' @name coerce,SingleCellExperiment,seurat-method
#'
#' @section `SingleCellExperiment` to `seurat`:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `seurat` class object. Use
#' [convertGenesToSymbols()] to convert gene IDs to names (symbols).
setAs(
    from = "SingleCellExperiment",
    to = "seurat",
    def = as.seurat.SingleCellExperiment
)



# seurat to SingleCellExperiment ===============================================
as.SingleCellExperiment <- function(from) {
    UseMethod("as.SingleCellExperiment")
}



# Fast internal method, that ensures dimnames are sanitized.
# Sanitizes legacy seurat objects that contain Ensembl IDs in names.
as.SingleCellExperiment.seurat <- function(from) {  # nolint
    validObject(from)
    names(rownames(from@raw.data)) <- NULL
    names(rownames(from@data)) <- NULL
    Seurat::Convert(from = from, to = "sce")
}



#' @rdname as
#' @name coerce,seurat,SingleCellExperiment-method
#'
#' @section `seurat` to `SingleCellExperiment`:
#' S4 coercion support for creating a `SingleCellExperiment` from a `seurat`
#' class object. The [Seurat FAQ page](https://satijalab.org/seurat/faq)
#' explains the `seurat` S4 class structure in detail. Internally, this method
#' improves the basic `Seurat::as.SingleCellExperiment()` S3 coercion method,
#' including the `object@scale.data` matrix, and will keep track of stashed
#' `rowRanges` and `metadata` if the `seurat` object was originally created
#' from a `SingleCellExperiment` (i.e. from the bcbioSingleCell package).
setAs(
    from = "seurat",
    to = "SingleCellExperiment",
    def = function(from) {
        to <- as.SingleCellExperiment(from)

        # Assays.
        assays <- assays(to)
        # Slot scaleData, if defined. Note that we're making the dimensions
        # match the other count matrices here.
        scaleData <- slot(from, "scale.data")
        if (!is.null(scaleData)) {
            scaleData <- scaleData[rownames(to), colnames(to)]
            assays[["scaleData"]] <- scaleData
        }
        stash <- .getSeuratStash(from, "assays")
        if (!is.null(stash)) {
            assert_are_disjoint_sets(names(assays), names(stash))
            assays <- c(assays, stash)
            assays(to) <- assays
        }

        # Row and column data.
        rowRanges(to) <- rowRanges(from)
        colData(to) <- colData(to)

        # Metadata.
        metadata(to) <- metadata(from)
        metadata(to)[["varGenes"]] <- slot(from, "var.genes")

        to
    }
)



# seurat to SummarizedExperiment ===============================================
#' @rdname as
#' @name coerce,seurat,RangedSummarizedExperiment-method
#'
#' @section `seurat` to `RangedSummarizedExperiment`:
#' S4 coercion support for creating a `RangedSummarizedExperiment` from a
#' `seurat` class object.
setAs(
    from = "seurat",
    to = "RangedSummarizedExperiment",
    def = function(from) {
        sce <- as(from, "SingleCellExperiment")
        rse <- as(sce, "RangedSummarizedExperiment")
        rse
    }
)



#' @rdname as
#' @name coerce,seurat,SummarizedExperiment-method
#' @section `seurat` to `SummarizedExperiment`:
#' S4 coercion support for creating a `SummarizedExperiment` from a `seurat`
#' class object.
setAs(
    from = "seurat",
    to = "SummarizedExperiment",
    def = function(from) {
        rse <- as(from, "RangedSummarizedExperiment")
        se <- as(rse, "SummarizedExperiment")
        se
    }
)



# Markers to tbl_df ======================================================
#' @rdname as
#' @name coerce,Markers,tbl_df-method
#' @section `Markers` to `tbl_df`:
#' S4 coercion support for creating a `tbl_df` from a `Markers` object.
setAs(
    from = "Markers",
    to = "tbl_df",
    def = function(from) {
        validObject(from)

        # Get gene2symbol from slotted ranges.
        g2s <- mcols(from[["ranges"]])[c("geneID", "geneName")]

        data <- as(from, "DataFrame")
        data[["ranges"]] <- NULL
        assert_are_disjoint_sets(colnames(data), colnames(g2s))
        data <- cbind(data, g2s)
        data <- as(data, "tbl_df")

        # Inform the user when sanitizing `Seurat::FindAllMarkers()` return.
        if ("cluster" %in% colnames(data)) {
            message("Grouping by cluster.")
            data <- group_by(data, !!sym("cluster"))
        }

        data
    }
)

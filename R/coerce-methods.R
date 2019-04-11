# FIXME Seurat v3 now masks these S3 methods.



#' Force an object to belong to a class
#'
#' @inherit methods::as
#' @name as
#' @aliases coerce
#' @importFrom methods coerce
#' @exportMethod coerce
#'
#' @seealso
#' - `Seurat::CreateSeuratObject()`.
#' - `Seurat::Convert()`.
#'
#' @examples
#' data(sce, package = "acidtest")
#' data(seurat_small)
#'
#' ## SingleCellExperiment to Seurat ====
#' object <- sce
#' print(object)
#' x <- as(object, "Seurat")
#' class(x)
#' print(x)
#'
#' ## Seurat to SingleCellExperiment ====
#' x <- as(seurat_small, "SingleCellExperiment")
#' print(x)
NULL



# CellCycleMarkers to tbl_df ===================================================
#' @rdname as
#' @name coerce,CellCycleMarkers,tbl_df-method
#' @section `CellCycleMarkers` to `tbl_df`:
#' S4 coercion support for creating a `tbl_df` from `CellCycleMarkers`.
setAs(
    from = "CellCycleMarkers",
    to = "tbl_df",
    def = function(from) {
        validObject(from)
        data <- do.call(what = rbind, args = from)
        data <- as(data, "tbl_df")
        message("Grouping by phase.")
        data <- group_by(data, !!sym("phase"))
        data
    }
)



# CellTypeMarkers to tbl_df ====================================================
#' @rdname as
#' @name coerce,CellTypeMarkers,tbl_df-method
#' @section `CellTypeMarkers` to `tbl_df`:
#' S4 coercion support for creating a `tbl_df` from `CellTypeMarkers`.
setAs(
    from = "CellTypeMarkers",
    to = "tbl_df",
    def = function(from) {
        validObject(from)
        data <- do.call(what = rbind, args = from)
        data <- as(data, "tbl_df")
        message("Grouping by cellType.")
        data <- group_by(data, !!sym("cellType"))
        data
    }
)



# SingleCellExperiment to Seurat ===============================================
# FIXME Seurat 3 now masks this method.
# Seurat:::as.Seurat.SingleCellExperiment
.as.Seurat.SingleCellExperiment <- function(from) {
    # Create the Seurat object.
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
    assert(identical(dim(from), dim(slot(to, "raw.data"))))

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



# FIXME Rethink this

#' @rdname as
#' @name coerce,SingleCellExperiment,Seurat-method
#'
#' @section `SingleCellExperiment` to `Seurat`:
#' Interally `Seurat::CreateSeuratObject` is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `Seurat` class object. Use
#' `convertGenesToSymbols` to convert gene IDs to names (symbols).
setAs(
    from = "SingleCellExperiment",
    to = "Seurat",
    def = .as.Seurat.SingleCellExperiment
)



# Seurat to SingleCellExperiment ===============================================
# Fast internal method, that ensures dimnames are sanitized.
# Sanitizes legacy Seurat objects that contain Ensembl IDs in names.
# FIXME This should be safe to deprecate.
# Seurat:::as.SingleCellExperiment.Seurat
.as.SingleCellExperiment.Seurat <- function(from) {  # nolint
    validObject(from)
    names(rownames(from@raw.data)) <- NULL
    names(rownames(from@data)) <- NULL
    Seurat::Convert(from = from, to = "sce")
}



# FIXME Rethink this approach

#' @rdname as
#' @name coerce,Seurat,SingleCellExperiment-method
#'
#' @section `Seurat` to `SingleCellExperiment`:
#' S4 coercion support for creating a `SingleCellExperiment` from a `Seurat`
#' class object. The [Seurat FAQ page](https://satijalab.org/Seurat/faq)
#' explains the `Seurat` S4 class structure in detail. Internally, this method
#' improves the basic `Seurat::as.SingleCellExperiment` S3 coercion method,
#' including the `object@scale.data` matrix, and will keep track of stashed
#' `rowRanges` and `metadata` if the `Seurat` object was originally created
#' from a `SingleCellExperiment` (i.e. from the bcbioSingleCell package).
setAs(
    from = "Seurat",
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
            assert(areDisjointSets(names(assays), names(stash)))
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



# Seurat to SummarizedExperiment ===============================================
#' @rdname as
#' @name coerce,Seurat,RangedSummarizedExperiment-method
#'
#' @section `Seurat` to `RangedSummarizedExperiment`:
#' S4 coercion support for creating a `RangedSummarizedExperiment` from a
#' `Seurat` class object.
setAs(
    from = "Seurat",
    to = "RangedSummarizedExperiment",
    def = function(from) {
        sce <- as(from, "SingleCellExperiment")
        rse <- as(sce, "RangedSummarizedExperiment")
        rse
    }
)



#' @rdname as
#' @name coerce,Seurat,SummarizedExperiment-method
#' @section `Seurat` to `SummarizedExperiment`:
#' S4 coercion support for creating a `SummarizedExperiment` from a `Seurat`
#' class object.
setAs(
    from = "Seurat",
    to = "SummarizedExperiment",
    def = function(from) {
        rse <- as(from, "RangedSummarizedExperiment")
        se <- as(rse, "SummarizedExperiment")
        se
    }
)



# SeuratMarkers to tbl_df ======================================================
#' @rdname as
#' @name coerce,SeuratMarkers,tbl_df-method
#' @section `SeuratMarkers` to `tbl_df`:
#' S4 coercion support for creating a `tbl_df` from a `Markers` object.
setAs(
    from = "SeuratMarkers",
    to = "tbl_df",
    def = function(from) {
        validObject(from)

        # Get gene2symbol from slotted ranges.
        g2s <- mcols(from[["ranges"]])[c("geneID", "geneName")]

        data <- as(from, "DataFrame")
        data[["ranges"]] <- NULL
        assert(areDisjointSets(colnames(data), colnames(g2s)))
        data <- cbind(data, g2s)
        data <- as(data, "tbl_df")

        data
    }
)



# SeuratMarkersPerCluster to tbl_df ============================================
#' @rdname as
#' @name coerce,SeuratMarkersPerCluster,tbl_df-method
#' @section `SeuratMarkersPerCluster` to `tbl_df`:
#' S4 coercion support for creating a `tbl_df` from `SeuratMarkersPerCluster`.
setAs(
    from = "SeuratMarkersPerCluster",
    to = "tbl_df",
    def = function(from) {
        validObject(from)
        data <- do.call(what = rbind, args = from)

        # Get gene2symbol from slotted ranges.
        g2s <- mcols(data[["ranges"]])[c("geneID", "geneName")]
        data[["ranges"]] <- NULL

        assert(areDisjointSets(colnames(data), colnames(g2s)))
        data <- cbind(data, g2s)
        # Ensure Rle columns get decoded.
        data <- decode(data)
        data <- as(data, "tbl_df")

        message("Grouping by cluster.")
        data <- group_by(data, !!sym("cluster"))

        data
    }
)

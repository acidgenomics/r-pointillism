#' Force an object to belong to a class
#'
#' @name coerce
#' @importFrom methods coerce
#' @exportMethod coerce
#'
#' @seealso
#' - `Seurat::CreateSeuratObject()`.
#' - `Seurat::as.Seurat()`.
#' - `Seurat::as.SingleCellExperiment()`.
#'
#' @examples
#' data(sce, package = "acidtest")
#' data(seurat)
#'
#' ## SingleCellExperiment to Seurat ====
#' object <- sce
#' print(object)
#' x <- as(object, "Seurat")
#' class(x)
#' print(x)
#'
#' ## Seurat to SingleCellExperiment ====
#' x <- as(seurat, "SingleCellExperiment")
#' print(x)
NULL



# CellCycleMarkers to tbl_df ===================================================
#' @rdname coerce
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
#' @rdname coerce
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



# Seurat to SingleCellExperiment ===============================================
#' @rdname coerce
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
        # Using the Seurat S3 coercion method here.
        to <- as.SingleCellExperiment(x = from, assay = NULL)

        # Row and column data.
        rowRanges(to) <- rowRanges(from)
        colData(to) <- colData(to)

        # Metadata.
        metadata(to) <- metadata(from)
        metadata(to)[["scaleData"]] <- GetAssayData(from, slot = "scale.data")
        metadata(to)[["variableFeatures"]] <- VariableFeatures(from)

        to
    }
)



# Seurat to RangedSummarizedExperiment =========================================
#' @rdname coerce
#' @name coerce,Seurat,RangedSummarizedExperiment-method
#'
#' @section `Seurat` to `RangedSummarizedExperiment`:
#' S4 coercion support for creating a `RangedSummarizedExperiment` from a
#' `Seurat` class object.
setAs(
    from = "Seurat",
    to = "RangedSummarizedExperiment",
    def = function(from) {
        to <- from
        to <- as(to, "SingleCellExperiment")
        to <- as(to, "RangedSummarizedExperiment")
        to
    }
)



# Seurat to SummarizedExperiment ===============================================
#' @rdname coerce
#' @name coerce,Seurat,SummarizedExperiment-method
#' @section `Seurat` to `SummarizedExperiment`:
#' S4 coercion support for creating a `SummarizedExperiment` from a `Seurat`
#' class object.
setAs(
    from = "Seurat",
    to = "SummarizedExperiment",
    def = function(from) {
        to <- from
        to <- as(to, "RangedSummarizedExperiment")
        to <- as(to, "SummarizedExperiment")
        to
    }
)



# SeuratMarkers to tbl_df ======================================================
#' @rdname coerce
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
#' @rdname coerce
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



# SingleCellExperiment to Seurat ===============================================
#' @rdname coerce
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
    def = function(from) {
        # Create the Seurat object. Note that `as.Seurat()` method requires
        # `logcounts` to be defined in `assays()`, so we're using
        # `CreateSeuratObject()` here instead.
        to <- CreateSeuratObject(
            counts = counts(from),
            project = "pointillism",
            assay = "RNA",
            # Already applied filtering cutoffs for cells and genes.
            min.cells = 0L,
            min.features = 0L,
            names.field = 1L,
            names.delim = "_",
            meta.data = as.data.frame(colData(from))
        )

        # Check that the dimensions match exactly.
        assert(identical(x = dim(from), y = dim(to)))

        # Stash rowRanges.
        rowRanges <- rowRanges(from)

        # Stash metadata.
        metadata <- metadata(from)
        # Update the session information.
        metadata[["sessionInfo"]] <- session_info()

        # Seurat v3 still recommends using `misc` slot.
        misc <- list(
            rowRanges = rowRanges,
            metadata = metadata
        )
        misc <- Filter(Negate(is.null), misc)
        slot(to, name = "misc") <- misc

        to
    }
)

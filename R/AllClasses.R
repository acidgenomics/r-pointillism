# SingleCellMarkers ============================================================
#' `SingleCellMarkers` Class
#'
#' Class containing essential elements generated during differential expression
#' analysis with either Seurat, edgeR, or DESeq2. This class is essentially a
#' `list` with validity checks to ensure that the slotted `DataFrame` and
#' `GRanges` correspond.
#'
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @slot data `grouped_df`.
#' @slot rowRanges `GRanges`.
setClass(
    Class = "SingleCellMarkers",
    slots = c(
        data = "tbl_df",
        rowRanges = "GRanges"
    ),
    prototype = list(
        rowRanges = GRanges()
    )
)

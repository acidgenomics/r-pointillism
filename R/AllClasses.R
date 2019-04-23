.prototypeMetadata <- list(
    version = packageVersion("pointillism"),
    date = Sys.Date()
)



# CellCycleMarkers =============================================================
#' Cell-cycle markers
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in [`metadata()`][S4Vectors::metadata].
#'
#' @family S4 classes
#' @export
#'
#' @return `CellTypeMarkers`
setClass(
    Class = "CellCycleMarkers",
    contains = "CompressedSplitDataFrameList",
    validity = function(object) {
        validate(
            identical(
                x = lapply(object[[1L]], class),
                y = list(
                    phase = "factor",
                    geneID = "character",
                    geneName = "character"
                )
            ),
            isSubset(
                x = c("version", "organism", "ensemblRelease", "date"),
                y = names(metadata(object))
            )
        )
    }
)



# CellTypeMarkers ==============================================================
#' Cell-type markers
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in [`metadata()`][S4Vectors::metadata].
#'
#' @family S4 classes
#' @export
#'
#' @return `CellTypeMarkers`
setClass(
    Class = "CellTypeMarkers",
    contains = "CompressedSplitDataFrameList",
    validity = function(object) {
        validate(
            identical(
                x = lapply(object[[1L]], class),
                y = list(
                    cellType = "factor",
                    geneID = "character",
                    geneName = "character"
                )
            ),
            isSubset(
                x = c("version", "organism", "ensemblRelease", "date"),
                y = names(metadata(object))
            )
        )
    }
)



# Known Markers ================================================================
#' Known markers
#'
#' Class containing known markers detected.
#'
#' @family S4 classes
#' @export
#'
#' @return Results are grouped by `cellType` column and arranged by adjusted
#'   *P* value (`padj`).
setClass(
    Class = "KnownMarkers",
    contains = "DataFrame",
    validity = function(object) {
        validate(
            isSubset(
                x = c(
                    "cellType",
                    "cluster",
                    "geneID",
                    "geneName",
                    "name",
                    "padj",
                    "pvalue"
                ),
                y = colnames(object)
            ),
            isSubset(
                x = c("alpha", "date", "version"),
                y = names(metadata(object))
            )
        )
    }
)



# Seurat Markers ===============================================================
#' Seurat markers
#'
#' Class containing essential elements for marker gene analysis.
#'
#' @family S4 classes
#' @export
#'
#' @return Results are arranged by adjusted *P* value (`padj`).
setClass(
    Class = "SeuratMarkers",
    contains = "DataFrame",
    validity = function(object) {
        validate(
            hasRownames(object),
            identical(
                x = colnames(object),
                y = c(
                    "avgLogFC",
                    "pct1",
                    "pct2",
                    "pvalue",
                    "padj",
                    "ranges"
                )
            )
        )
    }
)



# SeuratMarkersPerCluster ======================================================
#' Seurat markers per cluster
#'
#' Class containing essential elements for marker per cluster analysis.
#'
#' @family S4 classes
#' @export
#'
#' @return Results are split per `cluster` and arranged by adjusted *P* value
#'   (`padj`).
setClass(
    Class = "SeuratMarkersPerCluster",
    contains = "CompressedSplitDataFrameList",
    validity = function(object) {
        validate(
            all(grepl("^cluster", names(object))),
            identical(
                x = colnames(object)[[1L]],
                y = c(
                    "cluster",
                    "name",
                    "avgLogFC",
                    "pct1",
                    "pct2",
                    "pvalue",
                    "padj",
                    "ranges"
                )
            )
        )
    }
)

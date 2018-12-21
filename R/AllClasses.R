# package_version
setOldClass(Classes = class(packageVersion("base")))

# session_info
setOldClass(Classes = "session_info")



.prototypeMetadata <- list(
    version = packageVersion("pointillism"),
    date = Sys.Date()
)



# CellCycleMarkers =============================================================
#' Cell-Cycle Markers
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in `metadata()`.
#'
#' @family S4 classes
#' @export
#'
#' @seealso `CellCycleMarkers()`.
#'
#' @return `CellTypeMarkers`
setClass(
    Class = "CellCycleMarkers",
    contains = "CompressedSplitDataFrameList"
)
setValidity(
    Class = "CellCycleMarkers",
    method = function(object) {
        validate_that(
            identical(
                x = lapply(object[[1L]], class),
                y = list(
                    phase = "factor",
                    geneID = "character",
                    geneName = "character"
                )
            ),
            is_subset(
                x = c("version", "organism", "ensemblRelease", "date"),
                y = names(metadata(object))
            )
        )
    }
)



# CellTypeMarkers ==============================================================
#' Cell-Type Markers
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in `metadata()`.
#'
#' @family S4 classes
#' @export
#'
#' @seealso `CellTypeMarkers()`.
#'
#' @return `CellTypeMarkers`
setClass(
    Class = "CellTypeMarkers",
    contains = "CompressedSplitDataFrameList"
)
setValidity(
    Class = "CellTypeMarkers",
    method = function(object) {
        validate_that(
            identical(
                x = lapply(object[[1L]], class),
                y = list(
                    cellType = "factor",
                    geneID = "character",
                    geneName = "character"
                )
            ),
            is_subset(
                x = c("version", "organism", "ensemblRelease", "date"),
                y = names(metadata(object))
            )
        )
    }
)



# Known Markers ================================================================
#' Known Markers
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
    contains = "DataFrame"  # FIXME SplitDataFrameList
)
setValidity(
    Class = "KnownMarkers",
    method = function(object) {
        validate_that(
            is_subset(
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
            is_subset(
                x = c("alpha", "date", "version"),
                y = names(metadata(object))
            )
        )
    }
)



# Seurat Markers ===============================================================
#' Seurat Markers
#'
#' Class containing essential elements for marker gene analysis.
#'
#' @family S4 classes
#' @export
#'
#' @return Results are arranged by adjusted *P* value (`padj`).
setClass(
    Class = "SeuratMarkers",
    contains = "DataFrame"
)
setValidity(
    Class = "SeuratMarkers",
    method = function(object) {
        validate_that(
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
#' Seurat Markers per Cluster
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
    contains = "CompressedSplitDataFrameList"
)
setValidity(
    Class = "SeuratMarkersPerCluster",
    method = function(object) {
        validate_that(
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

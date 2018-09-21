# S3 classes ===================================================================
# package_version
setOldClass(Classes = class(packageVersion("base")))

# session_info
setOldClass(Classes = "session_info")



# CellCycleMarkers =============================================================
#' `CellCycleMarkers` Class
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in [metadata()].
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
setClass(
    Class = "CellCycleMarkers",
    contains = "DataFrame"
)



# CellTypeMarkers ==============================================================
#' `CellTypeMarkers` Class
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in [metadata()].
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @return `CellTypeMarkers`.
setClass(
    Class = "CellTypeMarkers",
    contains = "DataFrame"
)



# SeuratMarkers ================================================================
#' `SeuratMarkers` Class
#'
#' Class containing essential elements for Seurat marker analysis.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
setClass(
    Class = "SeuratMarkers",
    contains = "DataFrame"
)

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
#'
#' @seealso [CellCycleMarkers()].
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
#' @seealso [CellTypeMarkers()].
setClass(
    Class = "CellTypeMarkers",
    contains = "DataFrame"
)



# SeuratMarkers ================================================================
#' `SeuratMarkers` Class
#'
#' Class containing essential elements for Seurat marker gene analysis.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @seealso [SeuratMarkers()].
setClass(
    Class = "SeuratMarkers",
    contains = "DataFrame"
)



# KnownSeuratMarkers ===========================================================
#' `KnownSeuratMarkers` Class
#'
#' Class containing essential elements for analysis of known, cell-type specific
#' marker genes.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @seealso [knownMarkers()].
setClass(
    Class = "KnownSeuratMarkers",
    contains = "SeuratMarkers"
)

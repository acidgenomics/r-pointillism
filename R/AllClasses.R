# S3 classes ===================================================================
# package_version
setOldClass(Classes = class(packageVersion("base")))

# session_info
setOldClass(Classes = "session_info")

# grouped_df
# Note that `tbl` class is giving S4 inheritance issues here.
setOldClass(Classes = c("grouped_df", "tbl_df", "data.frame"))



# SeuratMarkers ================================================================
#' `SeuratMarkers` Class
#'
#' Class containing essential elements for Seurat marker analysis.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @slot data `grouped_df`. Sanitized Seurat markers data.
#' @slot rowRanges `GRanges`.
#' @slot organism `string`. Full Latin organism name.
#' @slot ensemblRelease `scalar integer`. Ensembl release version.
#' @slot version `package_version`.
#' @slot date `Date`. Date the object was saved.
#' @slot sessionInfo `session_info`. [sessioninfo::session_info()] return.
#'
#' @return `SeuratMarkers`.
setClass(
    Class = "SeuratMarkers",
    slots = c(
        data = "tbl_df",
        rowRanges = "GRanges",
        organism = "character",
        ensemblRelease = "integer",
        version = "package_version",
        date = "Date",
        sessionInfo = "session_info"
    ),
    prototype = prototype(
        version = packageVersion("pointillism"),
        date = Sys.Date(),
        sessionInfo = session_info()
    )
)



# FIXME Move `sanitizeSeuratMarkers()` here, rename, and deprecate.
# Describe `FindAllMarkers()`/`FindMarkers()` to sanitizeSeuratMarkers`.



setValidity(
    Class = "SeuratMarkers",
    method = function(object) {
        # .assertIsKnownMarkers(object)
        # requiredCols <- c(
        #     "cellType",  # bcbio
        #     "cluster",   # Seurat
        #     "geneID",    # bcbio
        #     "geneName",  # bcbio
        #     "avgLogFC",  # Seurat v2.1
        #     "padj"       # Seurat v2.1
        # )
        # assert_is_subset(requiredCols, colnames(object))
        #
        #
        #
        # .isSanitizedMarkers <- function(
        #     object,
        #     package = "Seurat"
        # ) {
        #     package <- match.arg(package)
        #
        #     # General checks -----------------------------------------------------------
        #     if (!is(object, "grouped_df")) {
        #         return(FALSE)
        #     } else if (
        #         is.null(attr(object, "vars")) ||
        #         attr(object, "vars") != "cluster"
        #     ) {
        #         return(FALSE)
        #     } else if (!"geneID" %in% colnames(object)) {
        #         return(FALSE)
        #     }
        #
        #     # Package-specific checks --------------------------------------------------
        #     if (package == "Seurat") {
        #         # Check for `Seurat::FindAllMarkers()` return.
        #         # These columns are output in an inconsistent format, so we'll sanitize
        #         # into lowerCamelCase.
        #         seuratBlacklist <- c(
        #             "avg_diff",   # Legacy, now "avg_logFC"
        #             "avg_logFC",  # Renamed in v2.1
        #             "gene",
        #             "p_val",      # We'll rename to pvalue, matching DESeq2
        #             "p_val_adj",  # New in v2.1, we'll rename to padj, matching DESeq2
        #             "pct.1",
        #             "pct.2"
        #         )
        #         if (any(seuratBlacklist %in% colnames(object))) {
        #             return(FALSE)
        #         } else {
        #             return(TRUE)
        #         }
        #     }
        # }
        #
        #
        #
        TRUE
    }
)



# SCDE =========================================================================
# `SCDE` Class
#
# Single cell differential expression results.
#
# TODO Add this class, which is returned from `diffExp()`.
# TODO Keep track of zinb weights, etc.
# Data provenance!!!
# #' @slot zeroWeights `string`. Package name and version used to perform
# #' @slot generator `string`. Package name and version used to generate the
# #'   markers.



# CellCycleMarkers =============================================================
#' `CellCycleMarkers` Class
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in [attributes()].
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @slot organism `string`. Full Latin organism name.
#' @slot ensemblRelease `scalar integer`. Ensembl release version.
#' @slot version `package_version`.
#' @slot date `Date`. Date the object was saved.
#'
#' @return `CellCycleMarkers`.
setClass(
    Class = "CellCycleMarkers",
    contains = "grouped_df",
    slots = list(
        organism = "character",
        ensemblRelease = "integer",
        version = "package_version",
        date = "Date"
    ),
    prototype = prototype(
        version = packageVersion("pointillism"),
        date = Sys.Date()
    )
)



#' `CellCycleMarkers` Generator
#'
#' Currently designed for internal use by the pointillism package.
#'
#' @family S4 Generators
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param data `grouped_df`. Grouped by `phase` column.
#'
#' @return `CellCycleMarkers`.
CellCycleMarkers <-  # nolint
    function(
        data,
        organism,
        ensemblRelease
    ) {
        assert_is_a_string(organism)
        assert_is_an_integer(ensemblRelease)
        new(
            Class = "CellCycleMarkers",
            data,
            organism = organism,
            ensemblRelease = ensemblRelease
        )
    }



setValidity(
    Class = "CellCycleMarkers",
    method = function(object) {
        # assert_is_all_of(object, "grouped_df")
        # assert_has_rows(object)
        # assert_are_identical(
        #     x = colnames(object),
        #     y = c("phase", "geneID", "geneName")
        # )
        # assert_are_identical(group_vars(object), "phase")
        # assert_is_subset(
        #     x = c("version", "organism", "ensemblRelease", "date"),
        #     y = names(attributes(object))
        # )
        TRUE
    }
)



# CellTypeMarkers ==============================================================
#' `CellTypeMarkers` Class
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in [attributes()].
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @slot organism `string`. Full Latin organism name.
#' @slot ensemblRelease `scalar integer`. Ensembl release version.
#' @slot version `package_version`.
#' @slot date `Date`. Date the object was saved.
#'
#' @return `CellTypeMarkers`.
setClass(
    Class = "CellTypeMarkers",
    contains = "grouped_df",
    slots = list(
        organism = "character",
        ensemblRelease = "integer",
        version = "package_version",
        date = "Date"
    ),
    prototype = prototype(
        version = packageVersion("pointillism"),
        date = Sys.Date()
    )
)



#' `CellTypeMarkers` Generator
#'
#' Currently designed for internal use by the pointillism package.
#'
#' @family S4 Generators
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param data `grouped_df`. Grouped by `phase` column.
#'
#' @return `CellTypeMarkers`.
CellTypeMarkers <-  # nolint
    function(
        data,
        organism,
        ensemblRelease
    ) {
        assert_is_a_string(organism)
        assert_is_an_integer(ensemblRelease)
        new(
            Class = "CellTypeMarkers",
            data,
            organism = organism,
            ensemblRelease = ensemblRelease
        )
    }



setValidity(
    Class = "CellTypeMarkers",
    method = function(object) {
        # assert_is_all_of(object, "grouped_df")
        # assert_has_rows(object)
        # assert_are_identical(
        #     x = colnames(object),
        #     y = c("cellType", "geneID", "geneName")
        # )
        # assert_are_identical(group_vars(object), "cellType")
        # assert_is_subset(
        #     x = c("version", "organism", "ensemblRelease", "date"),
        #     y = names(attributes(object))
        # )
        TRUE
    }
)

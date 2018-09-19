# TODO Consider defining classes for all markers and pairwise DE separately.
# TODO Rename `SingleCellMarkers` to `SeuratMarkers`?



# package_version ==============================================================
setOldClass(Classes = class(packageVersion("base")))



# tibble =======================================================================
setOldClass(Classes = class(tibble()))
setOldClass(Classes = c("grouped_df", class(tibble())))



# CellCycleMarkers =============================================================
#' `CellCycleMarkers` Class
#'
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @slot version `package_version`.
#' @slot organism `string`. Full Latin organism name.
#' @slot ensemblRelease `scalar integer`. Ensembl release version.
#' @slot date `Date`. Date the object was saved.
#'
setClass(
    Class = "CellCycleMarkers",
    contains = "grouped_df",
    slots = list(
        version = "package_version",
        organism = "character",
        ensemblRelease = "integer",
        date = "Date"
    ),
    prototype = prototype(
        version = packageVersion("pointillism"),
        date = Sys.Date()
    )
)

#' `CellCycleMarkers` Generator
#'
#' @inheritParams general
#' @param data `grouped_df`. Grouped by `phase` column.
#'
#' @export
CellCycleMarkers <-  # nolint
    function(
        data,
        organism,
        ensemblRelease
    ) {
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
        assert_is_all_of(object, "grouped_df")
        assert_has_rows(object)
        assert_are_identical(
            x = colnames(object),
            y = c("phase", "geneID", "geneName")
        )
        assert_are_identical(group_vars(object), "phase")
        assert_is_subset(
            x = c("version", "organism", "ensemblRelease", "date"),
            y = names(attributes(object))
        )
        TRUE
    }
)



# CellTypeMarkers ==============================================================
#' `CellTypeMarkers` Class
#'
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
setClass(
    Class = "CellTypeMarkers",
    contains = "grouped_df"
)

setValidity(
    Class = "CellTypeMarkers",
    method = function(object) {
        assert_is_all_of(object, "grouped_df")
        assert_has_rows(object)
        assert_are_identical(
            x = colnames(object),
            y = c("cellType", "geneID", "geneName")
        )
        assert_are_identical(group_vars(object), "cellType")
        assert_is_subset(
            x = c("version", "organism", "ensemblRelease", "date"),
            y = names(attributes(object))
        )
        TRUE
    }
)



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
#' @slot generator `string`. Package name and version used to generate the
#'   markers.
setClass(
    Class = "SingleCellMarkers",
    slots = c(
        data = "tbl_df",
        rowRanges = "GRanges",
        zeroWeights = "character",
        caller = "character",
        version = "package_version",
        organism = "character",
        ensemblRelease = "integer",
        date = "Date"
    ),
    prototype = prototype(
        zeroWeights = character(),
        version = packageVersion("pointillism"),
        date = Sys.Date()
    )
)

# FIXME Require information on the Ensembl release used.

setValidity(
    Class = "SingleCellMarkers",
    method = function(object) {
        .assertIsKnownMarkers(object)
        requiredCols <- c(
            "cellType",  # bcbio
            "cluster",   # Seurat
            "geneID",    # bcbio
            "geneName",  # bcbio
            "avgLogFC",  # Seurat v2.1
            "padj"       # Seurat v2.1
        )
        assert_is_subset(requiredCols, colnames(object))
        TRUE
    }
)



.isSanitizedMarkers <- function(
    object,
    package = "Seurat"
) {
    package <- match.arg(package)

    # General checks -----------------------------------------------------------
    if (!is(object, "grouped_df")) {
        return(FALSE)
    } else if (
        is.null(attr(object, "vars")) ||
        attr(object, "vars") != "cluster"
    ) {
        return(FALSE)
    } else if (!"geneID" %in% colnames(object)) {
        return(FALSE)
    }

    # Package-specific checks --------------------------------------------------
    if (package == "Seurat") {
        # Check for `Seurat::FindAllMarkers()` return.
        # These columns are output in an inconsistent format, so we'll sanitize
        # into lowerCamelCase.
        seuratBlacklist <- c(
            "avg_diff",   # Legacy, now "avg_logFC"
            "avg_logFC",  # Renamed in v2.1
            "gene",
            "p_val",      # We'll rename to pvalue, matching DESeq2
            "p_val_adj",  # New in v2.1, we'll rename to padj, matching DESeq2
            "pct.1",
            "pct.2"
        )
        if (any(seuratBlacklist %in% colnames(object))) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    }
}

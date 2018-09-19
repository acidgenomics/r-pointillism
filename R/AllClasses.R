# S3 classes ===================================================================
# package_version
setOldClass(Classes = class(packageVersion("base")))

# session_info
setOldClass(Classes = "session_info")

# tibble
# Note that `tbl_df` is properly exported in v1.4.99.9000.
# Including `tbl` class here is causing an S4 inheritance error.
# setOldClass(Classes = c("grouped_df", "tbl_df", "data.frame"))



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
    contains = "DataFrame"
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
#' @param data `DataFrame`. Grouped by `phase` column.
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
            as(data, "DataFrame"),
            metadata = list(
                version = packageVersion("pointillism"),
                date = Sys.Date(),
                organism = organism,
                ensemblRelease = ensemblRelease
            )
        )
    }



setValidity(
    Class = "CellCycleMarkers",
    method = function(object) {
        assert_is_all_of(object, "DataFrame")
        assert_has_rows(object)
        assert_are_identical(
            x = colnames(object),
            y = c("phase", "geneID", "geneName")
        )
        assert_is_factor(object[["phase"]])
        assert_is_subset(
            x = c("version", "organism", "ensemblRelease", "date"),
            y = names(metadata(object))
        )
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
#' @return `CellTypeMarkers`.
setClass(
    Class = "CellTypeMarkers",
    contains = "DataFrame"
)

# FIXME Ensure cellType is factor.



#' `CellTypeMarkers` Generator
#'
#' Currently designed for internal use by the pointillism package.
#'
#' @family S4 Generators
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param data `DataFrame`. Grouped by `phase` column.
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
            as(data, "DataFrame"),
            metadata = list(
                version = packageVersion("pointillism"),
                date = Sys.Date(),
                organism = organism,
                ensemblRelease = ensemblRelease
            )
        )
    }



setValidity(
    Class = "CellTypeMarkers",
    method = function(object) {
        assert_is_all_of(object, "DataFrame")
        assert_has_rows(object)
        assert_are_identical(
            x = colnames(object),
            y = c("cellType", "geneID", "geneName")
        )
        assert_is_factor(object[["cellType"]])
        assert_is_subset(
            x = c("version", "organism", "ensemblRelease", "date"),
            y = names(attributes(object))
        )
        TRUE
    }
)



# SeuratMarkers ================================================================
#' `SeuratMarkers` Class
#'
#' Class containing essential elements for Seurat marker analysis.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @slot data `DataFrame`. Sanitized Seurat markers data.
#' @slot GRanges `GRanges`.
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
        data = "DataFrame",
        GRanges = "GRanges",
        metadata = "list"
    )
)



#' Sanitize Seurat Markers
#'
#' This generator function is designed to take the original return from a Seurat
#' marker analysis and add corresponding gene annotations.
#'
#' @note [Seurat::FindAllMarkers()] maps the counts matrix rownames correctly in
#'   the `gene` column, whereas [Seurat::FindMarkers()] maps them correctly in
#'   the rownames of the returned marker `data.frame`.
#'
#' @family Marker Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param data `data.frame`. Unmodified [Seurat::FindAllMarkers()] or
#'   [Seurat::FindMarkers()] return.
#' @param GRanges `GRanges`. Gene annotations. Names must correspond to the
#'   rownames defined in `seurat@data`. The function will automatically subset
#'   and arrange this accordingly.
#'
#' @return `SeuratMarkers`. Results are arranged by adjusted P value.
#' @export
#'
#' @examples
#' object <- seurat_small
#'
#' # `FindAllMarkers()` return.
#' invisible(capture.output(
#'     all_markers <- Seurat::FindAllMarkers(object)
#' ))
#' all_sanitized <- SeuratMarkers(
#'     data = all_markers,
#'     GRanges = rowRanges(object)
#' )
#' glimpse(all_sanitized)
#'
#' # `FindMarkers()` return.
#' invisible(capture.output(
#'     ident_3_markers <- Seurat::FindMarkers(
#'         object = object,
#'         ident.1 = "1",
#'         ident.2 = NULL
#'     )
#' ))
#' ident_3_sanitized <- SeuratMarkers(
#'     data = ident_3_markers,
#'     GRanges = rowRanges(object)
#' )
#' glimpse(ident_3_sanitized)
SeuratMarkers <- function(data, GRanges) {
    assert_is_data.frame(data)
    assertHasRownames(data)
    assert_is_all_of(GRanges, "GRanges")
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = colnames(mcols(GRanges))
    )

    # Sanitize Seurat return ---------------------------------------------------
    # Coerce to tibble.
    data <- as(data, "tbl_df")
    # Standardize with camel case.
    data <- camel(data)

    # Map the Seurat matrix rownames to `rownames` column in tibble.
    if ("cluster" %in% colnames(data)) {
        message("`Seurat::FindAllMarkers()` return detected")
        all <- TRUE
        assert_is_subset("gene", colnames(data))
        data <- data %>%
            mutate(rowname = NULL) %>%
            rename(name = !!sym("gene"))
    } else {
        message("`Seurat::FindMarkers()` return detected")
        all <- FALSE
        data <- rename(data, name = !!sym("rowname"))
    }

    # Update legacy columns.
    if ("avgDiff" %in% colnames(data)) {
        message(paste(
            "Renaming legacy `avgDiff` column to `avgLogFC`",
            "(changed in Seurat v2.1)"
        ))
        data[["avgLogFC"]] <- data[["avgDiff"]]
        data[["avgDiff"]] <- NULL
    }

    # Rename P value columns to match DESeq2 conventions.
    if ("pVal" %in% colnames(data)) {
        data[["pvalue"]] <- data[["pVal"]]
        data[["pVal"]] <- NULL
    }
    if ("pValAdj" %in% colnames(data)) {
        data[["padj"]] <- data[["pValAdj"]]
        data[["pValAdj"]] <- NULL
    }

    # Ensure that required columns are present.
    requiredCols <- c(
        "name",
        "pct1",
        "pct2",
        "avgLogFC",     # Seurat v2.1.
        "padj",
        "pvalue"        # Renamed from `p_val`.
    )
    assert_is_subset(requiredCols, colnames(data))

    if (isTRUE(all)) {
        # `cluster` is only present in `FindAllMarkers() return`.
        data <- data %>%
            select(!!!syms(c("cluster", "name")), everything()) %>%
            group_by(!!sym("cluster")) %>%
            arrange(!!sym("padj"), .by_group = TRUE)
    } else {
        data <- data %>%
            select(!!sym("name"), everything()) %>%
            arrange(!!sym("padj"))
    }

    # GRanges ----------------------------------------------------------------
    # Require that all of the markers are defined in GRanges.
    names <- sort(unique(data[["name"]]))
    assert_is_subset(names, names(GRanges))
    GRanges <- GRanges[names]

    # Return -------------------------------------------------------------------
    new(
        Class = "SeuratMarkers",
        data = as(data, "DataFrame"),
        GRanges = GRanges,
        metadata = list(
            version = packageVersion("pointillism"),
            date = Sys.Date(),
            sessionInfo = session_info(include_base = TRUE)
        )
    )
}



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
        #     if (!is(object, "DataFrame")) {
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

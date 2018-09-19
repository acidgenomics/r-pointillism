# S3 classes ===================================================================
# package_version
setOldClass(Classes = class(packageVersion("base")))

# session_info
setOldClass(Classes = "session_info")

# tibble
# Note that `tbl_df` is properly exported in v1.4.99.9000.
# Including `tbl` class here is causing an S4 inheritance error.
# setOldClass(Classes = c("grouped_df", "tbl_df", "data.frame"))



# Internal generator for either cell cycle or cell type markers.
.cellMarkers <- function(
    file,
    gs,
    ws = 1L,
    organism,
    ensemblRelease,
    gene2symbol = NULL,
    class = c("CellCycleMarkers", "CellTypeMarkers")
) {
    class <- match.arg(class)
    if (class == "CellCycleMarkers") {
        group <- "phase"
    } else if (class == "CellTypeMarkers") {
        group <- "cellType"
    }

    assert_is_a_string(organism)
    assertIsAnImplicitInteger(ensemblRelease)
    assert_is_any_of(gene2symbol, c("gene2symbol", "NULL"))

    # Local markers from local file or Google Sheets.
    if (
        (!missing(file) && !missing(gs)) ||
        (missing(file) && missing(gs))
    ) {
        stop("Specify either `file` or `gs` (googlesheet)")
    } else if (!missing(file)) {
        # CSV file mode.
        assert_is_a_string(file)
        data <- import(file)
    } else if (!missing(gs)) {
        # Google Sheet mode.
        assert_is_a_string(gs)
        assert_is_scalar(ws)
        ss <- gs_key(gs)
        data <- gs_read(ss = ss, ws = ws)
    }

    # Get gene-to-symbol mappings from Ensembl.
    if (is.null(gene2symbol)) {
        gene2symbol <- makeGene2symbolFromEnsembl(
            organism = organism,
            release = ensemblRelease
        )
    }
    assert_is_all_of(gene2symbol, "gene2symbol")
    # We're coercing to tibble here to join back to the markers data.
    gene2symbol <- gene2symbol %>%
        as("DataFrame") %>%
        set_rownames(NULL) %>%
        as("tbl_df")

    # Coerce to tibble and sanitize.
    data <- data %>%
        as("tbl_df") %>%
        camel() %>%
        select(!!!syms(c(group, "geneID"))) %>%
        .[complete.cases(.), , drop = FALSE]

    # Warn user about markers that aren't present in the gene2symbol.
    # This is useful for informing about putative markers that aren't expressed.
    setdiff <- setdiff(data[["geneID"]], gene2symbol[["geneID"]])
    if (length(setdiff)) {
        warning(paste(
            "Markers missing from gene2symbol (not expressed):",
            printString(setdiff),
            sep = "\n"
        ))
    }

    intersect <- intersect(
        x = data[["geneID"]],
        y = gene2symbol[["geneID"]]
    )
    assert_is_non_empty(intersect)

    data <- data %>%
        filter(!!sym("geneID") %in% !!intersect) %>%
        mutate(!!sym(group) := as.factor(!!sym(group))) %>%
        unique() %>%
        left_join(gene2symbol, by = "geneID") %>%
        group_by(!!sym(group)) %>%
        arrange(!!sym("geneName"), .by_group = TRUE)

    new(
        Class = class,
        as(data, "DataFrame"),
        metadata = list(
            version = packageVersion("pointillism"),
            date = Sys.Date(),
            organism = organism,
            ensemblRelease = ensemblRelease
        )
    )
}



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
#' @inheritParams CellTypeMarkers
#'
#' @return `CellCycleMarkers`.
#'
#' @examples
#' # Google Sheets method.
#' x <- CellCycleMarkers(
#'     gs = "1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw",
#'     ws = "Homo_sapiens",
#'     organism = "Homo sapiens",
#'     ensemblRelease = 92L
#' )
CellCycleMarkers <-  # nolint
    function() {
        do.call(
            what = .cellMarkers,
            args = matchArgsToDoCall(
                args = list(class = "CellCycleMarkers")
            )
        )
    }
f <- formals(.cellMarkers)
f <- f[setdiff(names(f), "class")]
formals(CellCycleMarkers) <- f



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



#' Cell Markers
#'
#' Local spreadsheets (CSV, Excel) and Google Sheets are currently supported.
#'
#' Must contain `phase` and `geneID` columns.
#'
#' @family Marker Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams googlesheets::gs_read
#' @param file `string` or `missing`. Gene markers file path (CSV or Excel).
#' @param gs `string` or `missing`. Google Sheets `sheet_key` identifier, which
#'   can be located with the [googlesheets::gs_ls()] function.
#' @param gene2symbol `gene2symbol` or `NULL`. Use the stashed gene-to-symbol
#'   mappings from a `SingleCellExperiment` or `seurat` object.
#'
#' @seealso
#' - [googlesheets::gs_key()].
#' - [googlesheets::gs_read()].
#'
#' @return `CellTypeMarkers`.
#' @export
#'
#' @examples
#' # Google Sheets method.
#' x <- CellTypeMarkers(
#'     gs = "1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0",
#'     ws = "Homo_sapiens",
#'     organism = "Homo sapiens",
#'     ensemblRelease = 92L
#' )
CellTypeMarkers <-  # nolint
    function() {
        do.call(
            what = .cellMarkers,
            args = matchArgsToDoCall(
                args = list(class = "CellTypeMarkers")
            )
        )
    }
f <- formals(.cellMarkers)
f <- f[setdiff(names(f), "class")]
formals(CellTypeMarkers) <- f



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
            y = names(metadata(object))
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
#'   the ranges and arrange them alphabetically.
#'
#' @return `SeuratMarkers`. Results are arranged by adjusted P value, and
#'   grouped per cluster if applicable.
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

    # Sanitize Seurat markers --------------------------------------------------
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
        data <- slot(object, name = "data")
        # `FindAllMarkers()`
        if ("cluster" %in% colnames(data)) {
            clusters <- data[["cluster"]]
            assert_is_factor(clusters)
            clusters <- levels(clusters)
            # Ensure that we haven't already defined `ident`.
            assert_are_identical(
                x = clusters,
                y = as.character(seq(from = 0L, to = length(clusters) - 1L))
            )
        }

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

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
            release = as.integer(ensemblRelease)
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
#' Cell-Cycle Markers
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
            args = matchArgsToDoCall(args = list(class = "CellCycleMarkers"))
        )
    }
f <- formals(.cellMarkers)
f <- f[setdiff(names(f), "class")]
formals(CellCycleMarkers) <- f



# CellTypeMarkers ==============================================================
#' Cell-Type Markers
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
            args = matchArgsToDoCall(args = list(class = "CellTypeMarkers"))
        )
    }
f <- formals(.cellMarkers)
f <- f[setdiff(names(f), "class")]
formals(CellTypeMarkers) <- f



# SeuratMarkers ================================================================
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
SeuratMarkers <- function(
    data,
    GRanges,
    alpha = 0.05
) {
    assert_is_data.frame(data)
    assert_has_rows(data)
    assertHasRownames(data)
    assert_is_all_of(GRanges, "GRanges")
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = colnames(mcols(GRanges))
    )
    assert_is_a_number(alpha)
    assert_all_are_in_open_range(alpha, lower = 0L, upper = 1L)

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
        as(data, "DataFrame"),
        metadata = list(
            alpha = alpha,
            GRanges = GRanges,
            version = packageVersion("pointillism"),
            date = Sys.Date(),
            sessionInfo = session_info(include_base = TRUE)
        )
    )
}

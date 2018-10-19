#' @inherit CellCycleMarkers-class
#' @export
CellCycleMarkers <- function(
    object,
    organism = NULL,
    ensemblRelease = NULL
) {
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



#' @inherit CellTypeMarkers-class
#' @export
CellTypeMarkers <- function(
    object,
    organism = NULL,
    ensemblRelease = NULL
) {
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



#' Import Markers
#'
#' Currently designed for internal use by the pointillism package.
#'
#' @section Google Sheets:
#'
#' Google Sheets requires OAuth authentication.
#'
#' See also:
#' - [googlesheets::gs_key()].
#' - [googlesheets::gs_read()].
#'
#' @name importMarkers
#'
#' @inheritParams general
#' @inheritParams basejump::makeGRanges
#' @inheritParams googlesheets::gs_read
#'
#' @param file `string`. File path (CSV or Excel).
#' @param gs `string`. Google Sheets `sheet_key` identifier, which
#'   can be located with the [googlesheets::gs_ls()] function.
#' @param gene2symbol `Gene2Symbol`. Gene-to-symbol mappings.
#'
#' @return `CellCycleMarkers` or `CellTypeMarkers`.
#'
#' @examples
#' organism <- "Homo sapiens"
#' genomeBuild <- "GRCh38"
#' ensemblRelease <- 92L
#'
#' ## Google sheets mode currently must be interactive.
#' if (isTRUE(interactive())) {
#'     library(googlesheets)
#'     gs_ls()
#'
#'     ## Gene-to-symbol mappings.
#'     g2s <- makeGene2SymbolFromEnsembl("Homo sapiens")
#'
#'     ## Cell-cycle markers.
#'     x <- importCellCycleMarkersFromGoogle(
#'         gs = "1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw",
#'         ws = "Homo_sapiens",
#'         gene2symbol = g2s
#'     )
#'
#'     ## Cell-type markers.
#'     x <- importCellTypeMarkersFromGoogle(
#'         gs = "1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0",
#'         ws = "Homo_sapiens",
#'         gene2symbol = g2s
#'     )
#' }
NULL



# Internal generator for either cell cycle or cell type markers.
.importMarkers <- function(
    file,
    gs,
    ws = 1L,
    gene2symbol = NULL,
    organism,
    genomeBuild,
    ensemblRelease,
    class = c("CellCycleMarkers", "CellTypeMarkers")
) {
    class <- match.arg(class)
    assert_is_all_of(gene2symbol, "Gene2Symbol")
    assert_is_a_string(organism)
    assert_is_a_string(genomeBuild)
    assertIsAnImplicitInteger(ensemblRelease)

    if (class == "CellCycleMarkers") {
        group <- "phase"
    } else if (class == "CellTypeMarkers") {
        group <- "cellType"
    }

    # Local markers from local file or Google Sheets.
    if (
        (!missing(file) && !missing(gs)) ||
        (missing(file) && missing(gs))
    ) {
        stop("Specify either `file` or `gs` (googlesheet).")
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

    # We're coercing to Gene2Symbol here to join back to the markers data.
    gene2symbol <- as_tibble(gene2symbol, rownames = NULL)

    # Coerce to tibble and sanitize.
    data <- data %>%
        as_tibble() %>%
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
        left_join(
            y = as_tibble(gene2symbol),
            by = "geneID"
        ) %>%
        group_by(!!sym(group)) %>%
        arrange(!!sym("geneName"), .by_group = TRUE)

    fun <- get(
        x = class,
        package = "pointillism",
        inherits = FALSE
    )
}



#' @describeIn importMarkers Cell-cycle markers.
#' @export
importCellCycleMarkersFromGoogle <-  # nolint
    function() {
        do.call(
            what = .importMarkers,
            args = matchArgsToDoCall(args = list(class = "CellCycleMarkers"))
        )
    }
f <- formals(.importMarkers)
f <- f[setdiff(names(f), c("class", "file"))]
formals(importCellCycleMarkersFromGoogle) <- f



#' @describeIn importMarkers Cell-type markers.
#' @export
importCellTypeMarkersFromGoogle <-  # nolint
    function() {
        do.call(
            what = .importMarkers,
            args = matchArgsToDoCall(
                args = list(class = "CellTypeMarkers")
            )
        )
    }
f <- formals(.importMarkers)
f <- f[setdiff(names(f), c("class", "file"))]
formals(importCellTypeMarkersFromGoogle) <- f

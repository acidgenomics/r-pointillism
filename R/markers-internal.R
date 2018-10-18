# Internal generator for either cell cycle or cell type markers.
.getMarkers <- function(
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
        as_tibble(rownames = NULL)

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

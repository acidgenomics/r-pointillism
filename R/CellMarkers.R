.cellMarkers <- function(
    object,
    gene2symbol,
    class = c("CellCycleMarkers", "CellTypeMarkers")
) {
    assert(
        is(object, "DataFrame"),
        is(gene2symbol, "Gene2Symbol")
    )
    class <- match.arg(class)

    if (class == "CellCycleMarkers") {
        group <- "phase"
    } else if (class == "CellTypeMarkers") {
        group <- "cellType"
    }

    # Coerce to tibble and sanitize.
    data <- object %>%
        as_tibble(rownames = NULL) %>%
        camel() %>%
        select(!!!syms(c(group, "geneID"))) %>%
        .[complete.cases(.), , drop = FALSE]

    # Warn user about markers that aren't present in the gene2symbol.
    # This is useful for informing about putative markers that aren't expressed.
    setdiff <- setdiff(data[["geneID"]], gene2symbol[["geneID"]])
    if (length(setdiff)) {
        stop(paste(
            "Markers missing from gene2symbol:",
            printString(setdiff),
            sep = "\n"
        ))
    }
    intersect <- intersect(data[["geneID"]], gene2symbol[["geneID"]])
    assert(isNonEmpty(intersect))

    data <- data %>%
        filter(!!sym("geneID") %in% !!intersect) %>%
        mutate(!!sym(group) := as.factor(!!sym(group))) %>%
        unique() %>%
        left_join(
            y = as_tibble(gene2symbol, rownames = NULL),
            by = "geneID"
        ) %>%
        group_by(!!sym(group)) %>%
        arrange(!!sym("geneName"), .by_group = TRUE) %>%
        as("DataFrame")

    out <- data %>%
        split(f = .[[group]], drop = FALSE) %>%
        snake()
    metadata(out) <- metadata(gene2symbol)
    out
}



#' @rdname CellCycleMarkers-class
#' @inheritParams basejump::params
#' @export
CellCycleMarkers <-  # nolint
    function(object, gene2symbol) {
        class <- "CellCycleMarkers"
        data <- .cellMarkers(
            object = object,
            gene2symbol = gene2symbol,
            class = class
        )
        new(Class = class, data)
    }



#' @rdname CellTypeMarkers-class
#' @inheritParams basejump::params
#' @export
CellTypeMarkers <-  # nolint
    function(object, gene2symbol) {
        class <- "CellTypeMarkers"
        data <- .cellMarkers(
            object = object,
            gene2symbol = gene2symbol,
            class = class
        )
        new(Class = class, data)
    }

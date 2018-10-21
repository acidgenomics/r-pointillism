#' @inherit CellCycleMarkers-class
#' @inheritParams general
#' @export
CellCycleMarkers <- function(object, gene2symbol) {
    class <- "CellCycleMarkers"
    data <- .prepareMarkers(
        object = object,
        gene2symbol = gene2symbol,
        class = class
    )
    new(Class = class, data)
}



#' @inherit CellTypeMarkers-class
#' @inheritParams general
#' @export
CellTypeMarkers <- function(object, gene2symbol) {
    class <- "CellTypeMarkers"
    data <- .prepareMarkers(
        object = object,
        gene2symbol = gene2symbol,
        class = class
    )
    new(Class = class, data)
}



.prepareMarkers <- function(
    object,
    gene2symbol,
    class = c("CellCycleMarkers", "CellTypeMarkers")
) {
    assert_is_all_of(object, "DataFrame")
    assert_is_all_of(gene2symbol, "Gene2Symbol")
    class <- match.arg(class)

    if (class == "CellCycleMarkers") {
        group <- "phase"
    } else if (class == "CellTypeMarkers") {
        group <- "cellType"
    }

    # Coerce to tibble and sanitize.
    data <- object %>%
        as_tibble() %>%
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
    intersect <- intersect(
        x = data[["geneID"]],
        y = gene2symbol[["geneID"]]
    )
    assert_is_non_empty(intersect)

    out <- data %>%
        filter(!!sym("geneID") %in% !!intersect) %>%
        mutate(!!sym(group) := as.factor(!!sym(group))) %>%
        unique() %>%
        left_join(
            y = as_tibble(gene2symbol),
            by = "geneID"
        ) %>%
        group_by(!!sym(group)) %>%
        arrange(!!sym("geneName"), .by_group = TRUE) %>%
        as("DataFrame")

    xxx <- split(x = out, f = out[["phase"]], drop = FALSE)
    xxx <- snake(xxx)
    assertHasValidNames(xxx)


    # Split these into a DataFrameList.

    metadata(out) <- metadata(gene2symbol)
    out
}

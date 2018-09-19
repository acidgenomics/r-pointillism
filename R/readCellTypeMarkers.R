# TODO Consider renaming to simple `CellTypeMarkers()` generator.



#' Read Cell Type Markers File
#'
#' @name readCellTypeMarkers
#' @family Marker Functions
#' @author Michael Steinbaugh
#'
#' @param file `string`. Gene markers file path (CSV or Excel).
#' @param gene2symbol `data.frame`. Gene-to-symbol mappings. Columns must
#'   contain `geneID` and `geneName`.
#'
#' @return `CellTypeMarkers`.
#' @export
#'
#' @examples
#' # Homo sapiens
#' file <- system.file(
#'     file.path("extdata", "cell_type_markers.csv"),
#'     package = "pointillism"
#' )
#' gene2symbol <- makeGene2symbolFromEnsembl("Homo sapiens")
#' x <- readCellTypeMarkers(file, gene2symbol = gene2symbol)
#' print(x)
readCellTypeMarkers <- function(file, gene2symbol) {
    assert_is_a_string(file)
    assert_is_all_of(gene2symbol, "gene2symbol")

    # Read the markers file as a tibble.
    data <- import(file) %>%
        as("tbl_df") %>%
        camel() %>%
        select(!!!syms(c("cellType", "geneID"))) %>%
        .[complete.cases(.), , drop = FALSE]

    # Coerce to tibble.
    gene2symbol <- gene2symbol %>%
        as("DataFrame") %>%
        set_rownames(NULL) %>%
        as("tbl_df")

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
        left_join(gene2symbol, by = "geneID") %>%
        unique() %>%
        group_by(!!sym("cellType")) %>%
        arrange(!!sym("geneName"), .by_group = TRUE)

    new("CellTypeMarkers", data)
}

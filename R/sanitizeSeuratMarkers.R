#' Sanitize Seurat Markers
#'
#' Currently only Seurat markers are supported. The function should properly
#' sanitize `seurat` objects using either gene IDs or gene symbols in the
#' rownames automatically.
#'
#' @note [Seurat::FindAllMarkers()] maps the counts matrix rownames correctly in
#'   the `gene` column, whereas [Seurat::FindMarkers()] maps them correctly in
#'   the rownames of the returned marker `data.frame`.
#'
#' @family Marker Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param data `data.frame`. [Seurat::FindAllMarkers()] or
#'   [Seurat::FindMarkers()] return.
#' @param rowRanges `GRanges`. Gene annotations. Must contain `geneName` column
#'   that corresponds to the gene names in the `seurat` object used to generate
#'   the markers `data.frame`.
#'
#' @return `SingleCellMarkers`. Results are arranged by adjusted P value.
#' @export
#'
#' @examples
#' object <- seurat_small
#'
#' # `FindAllMarkers()` return.
#' invisible(capture.output(
#'     all_markers <- Seurat::FindAllMarkers(object)
#' ))
#' all_sanitized <- sanitizeSeuratMarkers(
#'     data = all_markers,
#'     rowRanges = rowRanges(object)
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
#' ident_3_sanitized <- sanitizeSeuratMarkers(
#'     data = ident_3_markers,
#'     rowRanges = rowRanges(object)
#' )
#' glimpse(ident_3_sanitized)
sanitizeSeuratMarkers <- function(data, rowRanges) {
    assert_is_data.frame(data)
    assertHasRownames(data)
    assert_is_all_of(rowRanges, "GRanges")
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = colnames(mcols(rowRanges))
    )

    # Sanitize Seurat markers data.frame =======================================
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

    markers <- unique(data[["name"]])

    # rowRanges ================================================================
    # Require that all of the rownames are defined in rowRanges.
    assert_is_subset(markers, names(rowRanges))
    # Subset and arrange the GRanges to match.
    rowRanges <- rowRanges[markers]

    # Return SingleCellMarkers =================================================
    new(
        Class = "SingleCellMarkers",
        data = data,
        rowRanges = rowRanges
    )
}

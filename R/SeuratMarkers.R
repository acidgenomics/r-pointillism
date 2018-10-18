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
#'
#' @inheritParams general
#' @param markers `data.frame`. Unmodified [Seurat::FindAllMarkers()] or
#'   [Seurat::FindMarkers()] return.
#' @param ranges `GRanges`. Gene annotations. Names must correspond to the
#'   rownames defined in `seurat@data`. The function will automatically subset
#'   the ranges and arrange them alphabetically.
#'
#' @return `SeuratMarkers`. Results are arranged by adjusted P value, and
#'   grouped per cluster if applicable.
#' @export
#'
#' @examples
#' data(seurat_small)
#' object <- seurat_small
#'
#' # `FindAllMarkers()` return.
#' invisible(capture.output(
#'     all_markers <- Seurat::FindAllMarkers(object)
#' ))
#' all_sanitized <- SeuratMarkers(
#'     markers = all_markers,
#'     ranges = rowRanges(object)
#' )
#' str(all_sanitized)
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
#'     markers = ident_3_markers,
#'     ranges = rowRanges(object)
#' )
#' str(ident_3_sanitized)
SeuratMarkers <- function(
    markers,
    ranges,
    alpha = 0.05
) {
    assert_is_data.frame(markers)
    assert_has_rows(markers)
    assertHasRownames(markers)
    assert_is_all_of(ranges, "GRanges")
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = colnames(mcols(ranges))
    )
    assert_is_a_number(alpha)
    assert_all_are_in_open_range(alpha, lower = 0L, upper = 1L)

    # Sanitize Seurat markers --------------------------------------------------
    # Coerce to tibble.
    data <- as_tibble(markers, rownames = "rowname")
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

    # Bind ranges as column ---------------------------------------------------
    data <- as(data, "DataFrame")
    # Require that all of the markers are defined in ranges.
    assert_is_subset(unique(data[["name"]]), names(ranges))
    data[["ranges"]] <- ranges[data[["name"]]]

    # Return -------------------------------------------------------------------
    new(
        Class = "SeuratMarkers",
        data,
        metadata = list(
            alpha = alpha,
            version = packageVersion("pointillism"),
            date = Sys.Date(),
            sessionInfo = session_info(include_base = TRUE)
        )
    )
}

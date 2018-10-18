#' Cell Types per Cluster
#'
#' @name cellTypesPerCluster
#' @family Cluster Statistics Functions
#'
#' @inheritParams general
#' @param min `scalar integer`. Minimum number of marker genes per cluster.
#' @param max `scalar integer`. Maximum number of marker genes per cluster.
#'
#' @return `grouped_df`. Grouped by `cluster` column, containing the count (`n`)
#'   of significant known makers per cell type.
#'
#' @examples
#' data(known_markers_small)
#' x <- cellTypesPerCluster(known_markers_small)
#' print(x)
NULL



.cellTypesPerCluster.KnownSeuratMarkers <- function(
    object,
    min = 1L,
    max = Inf
) {
    validObject(object)
    assert_is_a_number(min)
    assert_is_a_number(max)
    assert_all_are_positive(c(min, max))

    # Note that the order is important here.
    groupCols <- c("cluster", "cellType")

    data <- as(object, "DataFrame")
    data[["geneID"]] <- mcols(data[["ranges"]])[["geneID"]]
    data[["geneName"]] <- mcols(data[["ranges"]])[["geneName"]]
    data[["ranges"]] <- NULL

    data <- data %>%
        as_tibble() %>%
        ungroup() %>%
        select(!!!syms(groupCols), everything()) %>%
        mutate_at(groupCols, as.factor) %>%
        group_by(!!!syms(groupCols)) %>%
        arrange(!!sym("padj"), .by_group = TRUE) %>%
        # Only positive markers are informative and should be used.
        filter(!!sym("avgLogFC") > 0L) %>%
        # Use `toString()` instead of `aggregate()` for R Markdown tables.
        # Genes are arranged by P value.
        summarize(
            n = n(),
            name = toString(!!sym("name")),
            geneID = toString(!!sym("geneID")),
            geneName = toString(!!sym("geneName"))
        ) %>%
        group_by(!!sym("cluster")) %>%
        arrange(desc(!!sym("n")), .by_group = TRUE)

    # Apply minimum and maximum gene cutoffs.
    if (is.numeric(min) && min > 1L) {
        data <- filter(data, !!sym("n") >= !!min)
    }
    if (is.numeric(max) && max > 1L) {
        data <- filter(data, !!sym("n") <= !!max)
    }
    assert_has_rows(data)

    data
}



#' @rdname cellTypesPerCluster
#' @export
setMethod(
    f = "cellTypesPerCluster",
    signature = signature("KnownSeuratMarkers"),
    definition = .cellTypesPerCluster.KnownSeuratMarkers
)

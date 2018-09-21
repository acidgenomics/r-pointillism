# FIXME Need to switch to S4 method.



#' Cell Types per Cluster
#'
#' @name cellTypesPerCluster
#' @family Cluster Statistics Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param data `data.frame` grouped by `cellType`. [knownMarkersDetected()]
#'   return.
#' @param min `scalar integer`. Minimum number of marker genes per cluster.
#' @param max `scalar integer`. Maximum number of marker genes per cluster.
#'
#' @return `grouped_df` grouped by "`cluster`", containing the count (`n`) of
#'   significant known makers per cell type.
#' @export
#'
#' @examples
#' x <- cellTypesPerCluster(known_markers_detected_small)
#' glimpse(x)
NULL



.cellTypesPerCluster <- function(
    data,
    min = 1L,
    max = Inf
) {
    .assertIsKnownMarkersDetected(data)
    assert_is_a_number(min)
    assert_is_a_number(max)

    # Note that the order is important here
    groupCols <- c("cluster", "cellType")

    data <- data %>%
        ungroup() %>%
        select(!!!syms(groupCols), everything()) %>%
        mutate_at(groupCols, as.factor) %>%
        group_by(!!!syms(groupCols)) %>%
        arrange(!!sym("padj"), .by_group = TRUE) %>%
        # Only positive markers are informative and should be used.
        filter(!!sym("avgLogFC") > 0L) %>%
        # Use `toString()` instead of `aggregate()` for R Markdown tables.
        summarize(
            n = n(),
            # Genes are arranged by P value
            geneID = toString(!!sym("geneID")),
            geneName = toString(!!sym("geneName")),
            rowname = toString(!!sym("rowname"))
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

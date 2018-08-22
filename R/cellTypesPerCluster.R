#' Cell Types per Cluster
#'
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param min `scalar integer`. Minimum number of marker genes per cluster.
#' @param max `scalar integer`. Maximum number of marker genes per cluster.
#'
#' @return `grouped_df` grouped by "`cluster`", containing the count (`n`) of
#'   significant known makers per cell type.
#' @export
#'
#' @examples
#' x <- cellTypesPerCluster(known_markers_small)
#' glimpse(x)
cellTypesPerCluster <- function(
    data,
    min = 1L,
    max = Inf
) {
    stopifnot(is(data, "grouped_df"))
    assert_are_identical(attr(data, "vars"), "cellType")
    requiredCols <- c(
        "cellType",  # bcbio
        "cluster",   # Seurat
        "geneID",    # bcbio
        "geneName",  # bcbio
        "avgLogFC",  # Seurat v2.1
        "padj"       # Seurat v2.1
    )
    assert_is_subset(requiredCols, colnames(data))

    # Note that the order is important here
    groupCols <- c("cluster", "cellType")

    data <- data %>%
        ungroup() %>%
        select(!!!syms(groupCols), everything()) %>%
        mutate_at(groupCols, as.factor) %>%
        group_by(!!!syms(groupCols)) %>%
        arrange(!!sym("padj"), .by_group = TRUE) %>%
        # Use `toString()` instead of `aggregate()` for R Markdown tables
        summarize(
            n = n(),
            # Genes are arranged by P value
            geneID = toString(!!sym("geneID")),
            geneName = toString(!!sym("geneName")),
            rowname = toString(!!sym("rowname"))
        ) %>%
        group_by(!!sym("cluster")) %>%
        arrange(desc(!!sym("n")), .by_group = TRUE)

    # Apply minimum and maximum gene cutoffs
    if (is.numeric(min) && min > 1L) {
        data <- filter(data, !!sym("n") >= !!min)
    }
    if (is.numeric(max) && max > 1L) {
        data <- filter(data, !!sym("n") <= !!max)
    }
    assert_has_rows(data)

    data
}

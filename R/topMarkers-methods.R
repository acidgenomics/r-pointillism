#' Top Markers
#'
#' @name topMarkers
#' @family Marker Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param n `scalar integer`. Number of genes per cluster.
#' @param direction `string`. Whether to include "`positive`", "`negative`", or
#'   "`both`" directions of association per cluster.
#' @param coding `boolean`. Include only protein coding genes.
#'
#' @seealso
#' - [dplyr::slice()].
#' - [dplyr::top_n()].
#'
#' @return `grouped_df`.
#'
#' @examples
#' data(all_markers_small)
#' x <- topMarkers(
#'     all_markers_small,
#'     n = 2L,
#'     direction = "positive"
#' )
#' print(x)
NULL



.topMarkers.SeuratMarkers <-  # nolint
    function(
        object,
        n = 10L,
        direction = c("positive", "negative", "both")
    ) {
        validObject(object)
        assertIsAnImplicitInteger(n)
        direction <- match.arg(direction)

        # FIXME Consider making this a shared internal function.
        data <- as(object, "DataFrame")
        data[["ranges"]] <- NULL
        data <- as(data, "tbl_df")

        if ("cluster" %in% colnames(data)) {
            message("Grouping by cluster")
            data <- group_by(data, !!sym("cluster"))
        }

        # Subset to positive or negative correlation, if desired ("direction")
        # Note that `avgDiff` has been renamed to `avgLogFC` in Seurat v2.1.
        if (direction == "positive") {
            data <- filter(data, !!sym("avgLogFC") > 0L)
        } else if (direction == "negative") {
            data <- filter(data, !!sym("avgLogFC") < 0L)
        }

        data %>%
            # Arrange by adjusted P value.
            arrange(!!sym("padj"), .by_group = TRUE) %>%
            # Take the top rows by using slice.
            dplyr::slice(1L:n)
    }



#' @rdname topMarkers
#' @export
setMethod(
    f = "topMarkers",
    signature = signature("SeuratMarkers"),
    definition = .topMarkers.SeuratMarkers
)

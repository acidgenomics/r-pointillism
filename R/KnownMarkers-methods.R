# FIXME
# @examples
# # FIXME Need to update this.
# \dontrun{
# data(all_markers_small, known_markers_small)
# x <- KnownMarkers(
#     markers = all_markers_small,
#     known = known_markers_small
# )
# head(x)
# }



#' @inherit KnownMarkers-class
#' @name KnownMarkers
#'
#' @note Both the `markers` and `known` objects must contain Ensembl gene
#'   identifiers in the `geneID` column. We must avoid any matching operations
#'   based on the gene names, since these change often and can mismatch
#'   easily.
#'
#' @inheritParams general
#' @param markers `SeuratMarkers` or `SeuratMarkersPerCluster`.
#' @param known `CellTypeMarkers`. Grouped by `cellType` column. Known markers
#'   `data.frame` imported by [readCellTypeMarkers()] or pulled from internal
#'   cell cycle markers data.
#' @param promiscuousThreshold `scalar integer`. Minimum number of clusters
#'   required to consider a gene marker promiscuous. Set to `0` to disable
#'   promiscuous marker filtering.
#'
#' @return `KnownMarkers`.
NULL



KnownMarkers.SeuratMarkersPerCluster <-  # nolint
    function(
        markers,
        known,
        promiscuousThreshold = 0L
    ) {
        validObject(markers)
        validObject(known)
        assertIsAnImplicitInteger(promiscuousThreshold)
        assert_all_are_non_negative(promiscuousThreshold)
        promiscuousThreshold <- as.integer(promiscuousThreshold)

        alpha <- metadata(markers)[["alpha"]]
        assertIsAlpha(alpha)

        # Coerce markers data to tibble.
        markers <- do.call(rbind, markers)
        known <- do.call(rbind, known)

        # Determine where the known markers are located in the markers data.
        # Here we have slotted the gene IDs inside a "ranges" column.
        allGenes <- unique(data[["geneID"]])
        knownGenes <- known[["geneID"]]
        assert_are_intersecting_sets(knownGenes, allGenes)
        keep <- allGenes %in% knownGenes
        data <- markers[keep, , drop = FALSE]

        # Apply our alpha level cutoff.
        keep <- data[["padj"]] < alpha
        data <- data[keep, , drop = FALSE]

        # Add the `cellType` column.
        map <- left_join(
            x = tibble(
                name = data[["name"]],
                geneID = mcols(data[["ranges"]])[["geneID"]]
            ),
            y = known %>%
                as_tibble(rownames = NULL) %>%
                select(!!!syms(c("cellType", "geneID"))),
            by = "geneID"
        )
        assert_are_identical(data[["name"]], map[["name"]])
        data[["cellType"]] <- map[["cellType"]]

        # Filter out promiscuous markers present in multiple clusters.
        if (promiscuousThreshold > 1L) {
            cols <- c("cellType", "name")
            promiscuous <- data[, cols] %>%
                as_tibble() %>%
                ungroup() %>%
                group_by(!!!syms(cols)) %>%
                summarize(n = n()) %>%
                filter(!!sym("n") >= !!promiscuousThreshold) %>%
                pull("name")
            if (length(promiscuous)) {
                message(paste(
                    "Removing promiscuous markers:", toString(promiscuous)
                ))
                keep <- !data[["name"]] %in% promiscuous
                data <- data[keep, , drop = FALSE]
            }
        }

        metadata(data) <- list(
            alpha = alpha,
            version = packageVersion("pointillism"),
            date = Sys.Date()
        )
        new(Class = "KnownMarkers", data)
    }



#' @rdname KnownMarkers
#' @export
setMethod(
    f = "KnownMarkers",
    signature = signature(
        markers = "SeuratMarkersPerCluster",
        known = "CellTypeMarkers"
    ),
    definition = KnownMarkers.SeuratMarkersPerCluster
)

# FIXME Rethink how we're approaching this.



#' Known Markers Detected
#'
#' @note Both the `all` and `known` objects must contain Ensembl gene
#'   identifiers in the `geneID` column. We must avoid any matching operations
#'   based on the gene names, since these change often and can mismatch
#'   easily.
#'
#' @name knownMarkers
#' @family Marker Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param all `SeuratMarkers`.
#' @param known `CellTypeMarkers`. Grouped by `cellType` column. Known markers
#'   `data.frame` imported by [readCellTypeMarkers()] or pulled from internal
#'   [cell_cycle_markers] data.
#' @param promiscuousThreshold `scalar integer`. Minimum number of clusters
#'   required to consider a gene marker promiscuous. Set to `0` to disable
#'   promiscuous marker filtering.
#'
#' @return `SeuratMarkers`.
#'
#' @examples
#' x <- knownMarkers(
#'     all = all_markers_small,
#'     known = known_markers_small
#' )
#' head(x)
NULL



.knownMarkers.SeuratMarkers <- function(
    all,
    known,
    promiscuousThreshold = 5L
) {
    assert_is_all_of(all, "SeuratMarkers")
    assert_is_all_of(known, "CellTypeMarkers")
    assertIsAnImplicitInteger(promiscuousThreshold)
    assert_all_are_non_negative(promiscuousThreshold)
    promiscuousThreshold <- as.integer(promiscuousThreshold)

    alpha <- metadata(all)[["alpha"]]
    assertIsAlpha(alpha)

    # Determine where the known markers are located in the all markers data.
    # Here we have slotted the gene IDs inside a "ranges" column.
    allGenes <- all %>%
        .[["ranges"]] %>%
        mcols() %>%
        .[["geneID"]]
    knownGenes <- known[["geneID"]]
    assert_are_intersecting_sets(knownGenes, allGenes)
    keep <- allGenes %in% knownGenes
    data <- all[keep, , drop = FALSE]

    # Apply our alpha level cutoff.
    keep <- data[["padj"]] < alpha
    data <- data[keep, , drop = FALSE]

    # Now add the cell type column.
    # How do we want to map and group the cell type?
    # FIXME
    stop("This is still in progress.")

    x <- DataFrame(
        order = seq_len(nrow(data)),
        name = data[["name"]],
        geneID = mcols(data[["ranges"]])[["geneID"]]
    )
    y <- as(known, "DataFrame")[, c("cellType", "geneID")]
    map <- merge(
        x = x,
        y = y,
        all.x = TRUE,
        by = "geneID"
    ) %>%
        .[.[["order"]], , drop = FALSE]
    data[["cellType"]] <- map[["cellType"]]

    # FIXME
    stop("Still a work in progress.")

    if (isTRUE(filterPromiscuous)) {
        # Filter out promiscuous markers present in multiple clusters.
        promiscuous <- data %>%
            ungroup() %>%
            group_by(!!!syms(c("cellType", "geneID", "geneName"))) %>%
            summarize(n = n()) %>%
            filter(!!sym("n") >= !!promiscuousThreshold) %>%
            pull("geneID")
        if (length(promiscuous)) {
            message(paste(
                "Promiscuous markers:", toString(promiscuous)
            ))
            data <- filter(data, !!sym("geneID") %in% !!!promiscuous)
        }
    }

    data
}

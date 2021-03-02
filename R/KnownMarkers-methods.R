#' Known markers
#'
#' @name KnownMarkers
#' @note Both the `markers` and `known` objects must contain Ensembl gene
#'   identifiers in the `geneId` column. We must avoid any matching operations
#'   based on the gene names, since these change often and can mismatch
#'   easily.
#' @note Updated 2020-10-12.
#'
#' @inheritParams AcidRoxygen::params
#' @param markers `SeuratMarkers` or `SeuratMarkersPerCluster`.
#' @param known `CellTypeMarkers`.
#'   Grouped by `cellType` column. Known markers `data.frame` imported by
#'   `readCellTypeMarkers` or pulled from internal cell cycle markers data.
#' @param promiscuousThreshold `integer(1)`.
#'   Minimum number of clusters required to consider a gene marker promiscuous.
#'   Set to `0` to disable promiscuous marker filtering.
#'
#' @return `KnownMarkers`.
#'
#' @examples
#' data(cell_type_markers_list, seurat_all_markers)
#'
#' ## SeuratMarkersPerCluster ====
#' markers <- seurat_all_markers
#' known <- cell_type_markers_list[["homoSapiens"]]
#' x <- KnownMarkers(markers = markers, known = known)
#' summary(x)
NULL



## Updated 2019-09-01.
`KnownMarkers,SeuratMarkersPerCluster,CellTypeMarkers` <-  # nolint
    function(
        markers,
        known,
        promiscuousThreshold = 0L
    ) {
        validObject(markers)
        validObject(known)
        assert(
            isInt(promiscuousThreshold),
            allAreNonNegative(promiscuousThreshold)
        )
        promiscuousThreshold <- as.integer(promiscuousThreshold)
        alphaThreshold <- metadata(markers)[["alphaThreshold"]]
        assert(isAlpha(alphaThreshold))
        markers <- unlist(markers, recursive = FALSE, use.names = FALSE)
        ranges <- markers[["ranges"]]
        markers[["ranges"]] <- NULL
        markers[["geneId"]] <- as.character(mcols(ranges)[["geneId"]])
        markers[["geneName"]] <- as.character(mcols(ranges)[["geneName"]])
        known <- unlist(known, recursive = FALSE, use.names = FALSE)
        known[["geneName"]] <- NULL
        ## Determine where the known markers are located in the markers data.
        ## Here we have slotted the gene IDs inside a "ranges" column.
        assert(areIntersectingSets(markers[["geneId"]], known[["geneId"]]))
        keep <- markers[["geneId"]] %in% known[["geneId"]]
        x <- markers[keep, , drop = FALSE]
        ## Apply our alpha level cutoff.
        keep <- x[["padj"]] < alphaThreshold
        x <- x[keep, , drop = FALSE]
        ## Add the `cellType` column.
        x <- leftJoin(x, known, by = "geneId")
        ## Filter out promiscuous markers present in multiple clusters.
        if (isTRUE(promiscuousThreshold > 1L)) {
            x <- .filterPromiscuousMarkers(x, n = promiscuousThreshold)
        }
        metadata(x) <- list(
            alphaThreshold = alphaThreshold,
            version = packageVersion(packageName()),
            date = Sys.Date()
        )
        new(Class = "KnownMarkers", x)
    }



#' @rdname KnownMarkers
#' @export
setMethod(
    f = "KnownMarkers",
    signature = signature(
        markers = "SeuratMarkersPerCluster",
        known = "CellTypeMarkers"
    ),
    definition = `KnownMarkers,SeuratMarkersPerCluster,CellTypeMarkers`
)

context("Plot Functions")



# plotCellTypesPerCluster ======================================================
test_that("plotCellTypesPerCluster : seurat", {
    cellTypesPerCluster <- cellTypesPerCluster(known_markers_small) %>%
        # Subset for speed
        head(2L)
    invisible(capture.output(
        p <- plotCellTypesPerCluster(
            object = sce_small,
            cellTypesPerCluster = cellTypesPerCluster
        )
    ))
    expect_is(p, "list")
    expect_is(p[[1L]][[1L]], "ggplot")
})



# plotFeature ==================================================================
test_that("plotFeatureTSNE : seurat", {
    p <- plotFeatureTSNE(sce_small, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})

test_that("plotFeatureUMAP : seurat", {
    p <- plotFeatureUMAP(sce_small, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})



# plotMarker ===================================================================
object <- sce_small
genes <- head(rownames(object))

test_that("plotMarkerTSNE : seurat", {
    args <- methodFormals("plotMarkerTSNE", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(args, function(arg) {
        p <- plotMarkerTSNE(object, genes = genes, expression = arg)
        expect_is(p, "ggplot")
    }))
})

test_that("plotMarkerUMAP : seurat", {
    args <- methodFormals("plotMarkerUMAP", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(args, function(arg) {
        p <- plotMarkerUMAP(object, genes = genes, expression = arg)
        expect_is(p, "ggplot")
    }))
})

test_that("plotKnownMarkersDetected : seurat", {
    invisible(capture.output(
        p <- plotKnownMarkersDetected(
            object = sce_small,
            markers = head(known_markers_small, 2L)
        )
    ))
    expect_is(p, "list")
})

test_that("plotTopMarkers : seurat", {
    markers <- topMarkers(all_markers_small, n = 1L) %>%
        # Subset for speed
        head(2L)
    invisible(capture.output(
        x <- plotTopMarkers(
            object = sce_small,
            markers = markers
        )
    ))
    expect_is(x, "list")
    expect_is(x[[1L]][[1L]], "ggplot")
})



# plotPCA ======================================================================
test_that("plotPCA : seurat", {
    p <- plotPCA(sce_small)
    expect_is(p, "ggplot")
})



# plotPCElbow ==================================================================
test_that("plotPCElbow : seurat", {
    x <- plotPCElbow(pbmc_small)
    expect_identical(x, seq_len(11L))
})



# plotTSNE =====================================================================
test_that("plotTSNE : seurat", {
    p <- plotTSNE(sce_small)
    expect_is(p, "ggplot")
})



# plotUMAP =====================================================================
test_that("plotUMAP : seurat", {
    p <- plotTSNE(sce_small)
    expect_is(p, "ggplot")
})

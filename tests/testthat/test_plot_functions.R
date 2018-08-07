context("Plot Functions")



# plotCellTypesPerCluster ======================================================
test_that("plotCellTypesPerCluster : seurat", {
    markers <- cellTypesPerCluster(known_markers_small) %>%
        # Subset for speed
        head(2L)
    invisible(capture.output(
        p <- plotCellTypesPerCluster(
            object = sce_small,
            markers = markers
        )
    ))
    expect_is(p, "list")
    expect_is(p[[1L]][[1L]], "ggplot")
})



# plotFeature ==================================================================
test_that("plotFeature : seurat", {
    p <- plotFeature(sce_small, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})



# plotMarker ===================================================================
object <- sce_small
genes <- head(rownames(object))

test_that("plotMarker : seurat", {
    expression <- methodFormals("plotMarker", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(expression, function(expression) {
        p <- plotMarker(object, genes = genes, expression = expression)
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
    x <- plotPCElbow(seurat_small)
    expect_identical(x, seq_len(17L))
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

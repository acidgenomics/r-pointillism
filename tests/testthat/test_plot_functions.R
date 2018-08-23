context("Plot Functions")



# plotCellTypesPerCluster ======================================================
test_that("plotCellTypesPerCluster", {
    markers <- known_markers_detected_small %>%
        cellTypesPerCluster() %>%
        # Subset for speed
        head(2L)
    invisible(capture.output(
        p <- plotCellTypesPerCluster(
            object = seurat_small,
            markers = markers
        )
    ))
    expect_is(p, "list")
    expect_is(p[[1L]][[1L]], "ggplot")
})



# plotFeature ==================================================================
test_that("plotFeature", {
    p <- plotFeature(sce_small, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})



# plotMarker ===================================================================
object <- sce_small
genes <- head(rownames(object))

test_that("plotMarker", {
    expression <- methodFormals("plotMarker", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(expression, function(expression) {
        p <- plotMarker(object, genes = genes, expression = expression)
        expect_is(p, "ggplot")
    }))
})

test_that("plotKnownMarkersDetected", {
    invisible(capture.output(
        p <- plotKnownMarkersDetected(
            object = sce_small,
            markers = head(known_markers_detected_small, 2L)
        )
    ))
    expect_is(p, "list")
})

test_that("plotTopMarkers", {
    markers <- topMarkers(all_markers_small, n = 1L) %>%
        # Subset for speed
        head(2L)
    invisible(capture.output(
        x <- plotTopMarkers(
            object = seurat_small,
            markers = markers
        )
    ))
    expect_is(x, "list")
    expect_is(x[[1L]][[1L]], "ggplot")
})



# plotPCA ======================================================================
test_that("plotPCA", {
    p <- plotPCA(sce_small)
    expect_is(p, "ggplot")
})



# plotPCElbow ==================================================================
test_that("plotPCElbow", {
    x <- plotPCElbow(Seurat::pbmc_small)
    expect_identical(x, seq_len(11L))
})



# plotTSNE =====================================================================
test_that("plotTSNE", {
    p <- plotTSNE(sce_small)
    expect_is(p, "ggplot")
})



# plotUMAP =====================================================================
test_that("plotUMAP", {
    p <- plotTSNE(sce_small)
    expect_is(p, "ggplot")
})

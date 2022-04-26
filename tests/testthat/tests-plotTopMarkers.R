context("plotTopMarkers")

test_that("Seurat", {
    object <- seurat
    markers <- seuratAllMarkers
    expect_s4_class(markers, "SeuratMarkersPerCluster")
    invisible(capture.output({
        x <- plotTopMarkers(
            object = object,
            markers = markers,
            direction = "up",
            n = 1L
        )
    }))
    expect_type(x, "list")
    expect_s3_class(x[[1L]][[1L]], "ggplot")
})

context("SeuratMarkers")

test_that("SeuratMarkers", {
    ranges <- rowRanges(seurat)
    invisible(capture.output(
        markers <- Seurat::FindMarkers(
            object = seurat,
            ident.1 = "1",
            ident.2 = NULL
        )
    ))
    x <- SeuratMarkers(object = markers, ranges = ranges)
    expect_s4_class(x, "SeuratMarkers")
})



context("SeuratMarkersPerCluster")

test_that("SeuratMarkersPerCluster", {
    ranges <- rowRanges(seurat)
    # Suppressing expected warning: "cannot compute exact p-value with ties"
    suppressWarnings(
        invisible(capture.output(
            markers <- Seurat::FindAllMarkers(seurat)
        ))
    )
    x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
    expect_s4_class(x, "SeuratMarkersPerCluster")
})

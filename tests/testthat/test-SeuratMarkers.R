context("SeuratMarkers")

test_that("SeuratMarkers", {
    object <- Seurat
    ranges <- rowRanges(object)
    invisible(capture.output({
        markers <- Seurat::FindMarkers(
            object = object,
            ident.1 = "1",
            ident.2 = NULL
        )
    }))
    x <- SeuratMarkers(object = markers, ranges = ranges)
    expect_s4_class(x, "SeuratMarkers")
})



context("SeuratMarkersPerCluster")

test_that("SeuratMarkersPerCluster", {
    object <- Seurat
    ranges <- rowRanges(object)
    ## Suppressing expected warning: "cannot compute exact p-value with ties".
    suppressWarnings({
        invisible(capture.output({
            markers <- Seurat::FindAllMarkers(object)
        }))
    })
    x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
    expect_s4_class(x, "SeuratMarkersPerCluster")
})

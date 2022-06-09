object <- objs[["Seurat"]]
ranges <- rowRanges(object)

test_that("SeuratMarkers", {
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



test_that("SeuratMarkersPerCluster", {
    ## Suppressing expected warning: "cannot compute exact p-value with ties".
    suppressWarnings({
        invisible(capture.output({
            markers <- Seurat::FindAllMarkers(object)
        }))
    })
    x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
    expect_s4_class(x, "SeuratMarkersPerCluster")
})

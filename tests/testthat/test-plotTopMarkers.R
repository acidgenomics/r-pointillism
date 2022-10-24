## Can hit this cryptic ggplot2 error on macOS unless we call at the top.
## > no applicable method for 'depth' applied to an object of class "NULL"
## https://github.com/tidyverse/ggplot2/issues/2514
grid::current.viewport()

test_that("Seurat", {
    object <- objs[["Seurat"]]
    markers <- objs[["SeuratMarkersPerCluster"]]
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

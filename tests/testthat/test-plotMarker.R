context("plotMarker")

expression <- eval(formals(`plotMarker,Seurat`)[["expression"]])
test_that("plotMarker", {
    for (object in objects) {
        invisible(lapply(
            X = expression,
            FUN = function(expression) {
                p <- plotMarker(
                    object = object,
                    genes = genes,
                    expression = expression
                )
                expect_s3_class(p, "ggplot")
            }
        ))
    }
})



context("plotKnownMarkers")

markers <- head(seuratKnownMarkers, n = 2L)
test_that("plotKnownMarkers", {
    for (object in objects) {
        invisible(capture.output({
            x <- plotKnownMarkers(
                object = object,
                markers = markers
            )
        }))
        expect_type(x, "list")
        expect_s3_class(x[[1L]][[1L]], "ggplot")
    }
})



context("plotTopMarkers")

test_that("Seurat", {
    object <- seurat
    markers <- head(seuratAllMarkers, n = 2L)
    invisible(capture.output({
        x <- plotTopMarkers(
            object = object,
            markers = markers
        )
    }))
    expect_type(x, "list")
    expect_s3_class(x[[1L]][[1L]], "ggplot")
})

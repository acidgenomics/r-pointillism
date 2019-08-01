context("plotMarker")

expression <- eval(formals(`plotMarker,Seurat`)[["expression"]])
with_parameters_test_that(
    "plotMarker", {
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
    },
    object = objects
)



context("plotKnownMarkers")

with_parameters_test_that(
    "plotKnownMarkers", {
        invisible(capture.output(
            p <- plotKnownMarkers(
                object = object,
                markers = seuratKnownMarkers
            )
        ))
        expect_type(p, "list")
    },
    object = objects
)



context("plotTopMarkers")

with_parameters_test_that(
    "plotTopMarkers", {
        markers <- head(seuratAllMarkers, n = 2L)
        invisible(capture.output(
            x <- plotTopMarkers(
                object = object,
                markers = markers
            )
        ))
        expect_type(x, "list")
        expect_s3_class(x[[1L]][[1L]], "ggplot")
    },
    object = objects
)

context("plotMarker")

with_parameters_test_that(
    "plotMarker", {
        expression <- methodFormals("plotMarker", "Seurat") %>%
            .[["expression"]] %>%
            as.character() %>%
            .[-1L]
        invisible(lapply(expression, function(expression) {
            p <- plotMarker(
                object = object,
                genes = genes,
                expression = expression
            )
            expect_s3_class(p, "ggplot")
        }))
    },
    object = objects
)



context("plotKnownMarkers")

with_parameters_test_that(
    "plotKnownMarkers", {
        invisible(capture.output(
            p <- plotKnownMarkers(
                object = object,
                markers = seurat_known_markers
            )
        ))
        expect_type(p, "list")
    },
    object = objects
)



context("plotTopMarkers")

with_parameters_test_that(
    "plotTopMarkers", {
        markers <- head(seurat_all_markers, n = 2L)
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

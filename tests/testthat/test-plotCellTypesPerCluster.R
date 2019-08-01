context("plotCellTypesPerCluster")

with_parameters_test_that(
    "plotCellTypesPerCluster", {
        invisible(capture.output(
            list <- plotCellTypesPerCluster(
                object = object,
                markers = seuratKnownMarkers
            )
        ))
        expect_type(list, "list")
        expect_s3_class(list[[1L]][[1L]], "ggplot")
    },
    object = objects
)

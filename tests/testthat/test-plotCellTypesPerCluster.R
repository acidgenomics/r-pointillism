context("plotCellTypesPerCluster")

## FIXME This unit test is now failing, following AcidExperiment update.
test_that("Homo sapiens", {
    for (object in objects) {
        invisible(capture.output({
            list <- plotCellTypesPerCluster(
                object = object,
                markers = seuratKnownMarkers
            )
        }))
        expect_type(list, "list")
        expect_s3_class(list[[1L]][[1L]], "ggplot")
    }
})

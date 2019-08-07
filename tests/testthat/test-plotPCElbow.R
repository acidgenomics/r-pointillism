context("plotPCElbow")

objects <- list(
    Seurat = seurat,
    cell_data_set = cds
)

with_parameters_test_that(
    "plotPCElbow", {
        x <- plotPCElbow(object)
        expect_identical(
            object = attr(x, "elbow"),
            expected = elbow
        )
    },
    object = objects,
    elbow = list(
        Seurat = 10L,
        cell_data_set = 16L
    )
)

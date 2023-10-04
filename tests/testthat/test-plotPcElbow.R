objects <- list(
    "Seurat" = objs[["Seurat"]]
    ## > "cell_data_set" = objs[["cell_data_set"]]
)
elbows <- list(
    "Seurat" = 10L
    ## > "cell_data_set" = 16L
)

test_that("plotPcElbow", {
    Map(
        object = objects,
        elbow = elbows,
        f = function(object, elbow) {
            x <- plotPcElbow(object)
            expect_identical(
                object = attr(x, "elbow"),
                expected = elbow
            )
        }
    )
})

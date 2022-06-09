objects <- list(
    "Seurat" = objs[["Seurat"]]
    ## > "cell_data_set" = objs[["cell_data_set"]]
)
elbows <- list(
    "Seurat" = 10L
    ## > "cell_data_set" = 16L
)

test_that("plotPCElbow", {
    Map(
        object = objects,
        elbow = elbows,
        f = function(object, elbow) {
            x <- plotPCElbow(object)
            expect_identical(
                object = attr(x, "elbow"),
                expected = elbow
            )
        }
    )
})

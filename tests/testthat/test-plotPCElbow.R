context("plotPCElbow")

objects <- list(
    "Seurat" = Seurat
    ## > "cell_data_set" = cell_data_set
)
elbows <- list(
    "Seurat" = 10L
    ## > "cell_data_set" = 16L
)

test_that("plotPCElbow", {
    mapply(
        object = objects,
        elbow = elbows,
        FUN = function(object, elbow) {
            x <- plotPCElbow(object)
            expect_identical(
                object = attr(x, "elbow"),
                expected = elbow
            )
        },
        SIMPLIFY = FALSE
    )
})

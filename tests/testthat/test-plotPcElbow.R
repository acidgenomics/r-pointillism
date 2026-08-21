objects <- list(
    Seurat = objs[["Seurat"]]
)
elbows <- list(
    Seurat = 10L
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

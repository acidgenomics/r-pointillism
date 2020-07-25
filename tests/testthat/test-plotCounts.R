context("plotCounts")

test_that("plotCounts", {
    for (object in objects) {
        ## Dot.
        p <- plotCounts(
            object = object,
            genes = genes,
            geom = "dot"
        )
        expect_s3_class(p, "ggplot")
        ## Violin.
        p <- plotCounts(
            object = object,
            genes = genes,
            geom = "violin"
        )
        expect_s3_class(p, "ggplot")
    }
})



context("plotDots")

test_that("plotDots", {
    for (object in objects) {
        p <- plotDots(object, genes = genes)
        expect_s3_class(p, "ggplot")
    }
})



context("plotViolin")

test_that("plotViolin", {
    for (object in objects) {
        p <- plotViolin(object, genes = genes)
        expect_s3_class(p, "ggplot")
    }
})

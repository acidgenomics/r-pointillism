context("plotReducedDim")

test_that("plotReducedDim", {
    for (object in objects) {
        p <- plotReducedDim(
            object = object,
            pointsAsNumbers = TRUE,
            dark = TRUE,
            label = FALSE
        )
        expect_s3_class(p, "ggplot")
    }
})

funs <- list(plotPCA, plotTSNE, plotUMAP)
test_that("Aliases", {
    for (fun in funs) {
        lapply(
            X = objects,
            FUN = function(object) {
                p <- fun(object)
                expect_s3_class(p, "ggplot")
            }
        )
    }
})

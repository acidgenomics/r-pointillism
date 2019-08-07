context("plotReducedDim")

with_parameters_test_that(
    "plotReducedDim", {
        p <- plotReducedDim(
            object = object,
            pointsAsNumbers = TRUE,
            dark = TRUE,
            label = FALSE
        )
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

funs <- list(plotPCA, plotTSNE, plotUMAP)

with_parameters_test_that(
    "Aliases", {
        lapply(
            X = objects,
            FUN = function(object) {
                p <- fun(object)
                expect_s3_class(p, "ggplot")
            }
        )
    },
    fun = funs
)

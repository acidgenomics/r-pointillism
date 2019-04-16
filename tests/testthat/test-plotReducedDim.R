context("plotReducedDim")

with_parameters_test_that(
    "plotReducedDim", {
        p <- plotReducedDim(
            object = object,
            reducedDim = "TSNE",
            pointsAsNumbers = TRUE,
            dark = TRUE,
            label = FALSE
        )
        expect_s3_class(p, "ggplot")
    },
    object = objects
)



context("plotPCA")

with_parameters_test_that(
    "plotPCA", {
        p <- plotPCA(object)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)



context("plotTSNE")

with_parameters_test_that(
    "plotTSNE", {
        p <- plotTSNE(object)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)



context("plotUMAP")

with_parameters_test_that(
    "plotUMAP", {
        p <- plotTSNE(object)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

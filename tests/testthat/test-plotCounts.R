context("plotCounts")

with_parameters_test_that(
    "plotCounts", {
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
    },
    object = objects
)



context("plotDot")

with_parameters_test_that(
    "plotDot", {
        p <- plotDot(object, genes = genes)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)



context("plotViolin")

with_parameters_test_that(
    "plotViolin", {
        p <- plotViolin(object, genes = genes)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

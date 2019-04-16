context("plotFeature")

with_parameters_test_that(
    "plotFeature", {
        features <- c("nCount_RNA", "nFeature_RNA")
        p <- plotFeature(object, features = features)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

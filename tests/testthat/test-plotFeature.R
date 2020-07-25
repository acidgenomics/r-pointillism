context("plotFeature")

test_that("plotFeature", {
    for (object in objects) {
        features <- c("nCount_RNA", "nFeature_RNA")
        p <- plotFeature(object, features = features)
        expect_s3_class(p, "ggplot")
    }
})

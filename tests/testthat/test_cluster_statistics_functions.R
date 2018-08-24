context("Cluster Statistics Functions")



# cellCountsPerCluster ========================================================
test_that("cellCountsPerCluster", {
    x <- cellCountsPerCluster(sce_small)
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "ident")
})



# cellTypesPerCluster ==========================================================
test_that("cellTypesPerCluster", {
    x <- cellTypesPerCluster(known_markers_detected_small)
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "cluster")
    expect_identical(
        lapply(x, class),
        list(
            cluster = "factor",
            cellType = "factor",
            n = "integer",
            geneID = "character",
            geneName = "character",
            rowname = "character"
        )
    )
})



# clusterCellCountsPerSample ==================================================
test_that("clusterCellCountsPerSample", {
    x <- clusterCellCountsPerSample(sce_small)
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "sampleName")
})

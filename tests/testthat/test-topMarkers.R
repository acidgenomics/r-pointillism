context("topMarkers")

test_that("Default", {
    object <- objs[["SeuratMarkersPerCluster"]]
    object <- topMarkers(object, direction = "up", n = 2L)
    expect_s4_class(object, "DataFrame")
    object <- lapply(X = object, FUN = class)
    object <- object[sort(names(object))]
    expect_identical(
        object = object,
        expected = list(
            "avgLog2Fc" = "numeric",
            "cluster" = "factor",
            "geneId" = "factor",
            "geneName" = "factor",
            "name" = "factor",
            "padj" = "numeric",
            "pct1" = "numeric",
            "pct2" = "numeric",
            "pvalue" = "numeric"
        )
    )
})

directions <- formals(`topMarkers,SeuratMarkersPerCluster`)[["direction"]]
test_that("direction", {
    for (direction in directions) {
        object <- objs[["SeuratMarkersPerCluster"]]
        object <- topMarkers(
            object = object,
            direction = direction
        )
        expect_s4_class(object, "DataFrame")
    }
})

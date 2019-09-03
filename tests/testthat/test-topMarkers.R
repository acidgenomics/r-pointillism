## FIXME Double check that this unit test is performing as expected.

context("topMarkers")

test_that("Default", {
    object <- topMarkers(seuratAllMarkers, direction = "up", n = 2L)
    expect_s4_class(object, "DataFrame")
    expect_identical(
        object %>% lapply(class) %>% .[sort(names(.))],
        list(
            avgLogFC = "numeric",
            cluster = "factor",
            geneID = "factor",
            geneName = "factor",
            name = "factor",
            padj = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            pvalue = "numeric"
        )
    )
})

direction <- formals(`topMarkers,SeuratMarkersPerCluster`)[["direction"]]
with_parameters_test_that(
    "direction", {
        object <- topMarkers(
            object = seuratAllMarkers,
            direction = direction
        )
        expect_s4_class(object, "DataFrame")
    },
    direction = direction
)

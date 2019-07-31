context("topMarkers")

test_that("grouped_df", {
    x <- topMarkers(seuratAllMarkers)
    expect_is(x, "grouped_df")
    expect_identical(group_vars(x), "cluster")
    expect_identical(
        lapply(x, class) %>% .[sort(names(.))],
        list(
            avgLogFC = "numeric",
            cluster = "factor",
            geneID = "character",
            geneName = "factor",
            name = "character",
            padj = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            pvalue = "numeric"
        )
    )
})

direction <- formals(topMarkers) %>%
    .[["direction"]] %>%
    as.character() %>%
    .[-1L]

with_parameters_test_that(
    "direction", {
        x <- topMarkers(
            data = seuratAllMarkers,
            direction = direction
        )
        expect_is(x, "tbl_df")
    },
    direction = direction
)

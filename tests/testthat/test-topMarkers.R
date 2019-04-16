context("topMarkers")

test_that("grouped_df", {
    x <- topMarkers(seurat_all_markers)
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "cluster")
    expect_identical(
        lapply(x, class) %>%
            .[sort(names(.))],
        list(
            avgLogFC = "numeric",
            cluster = "factor",
            geneID = "character",
            geneName = "character",
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

test_that("direction", {
    invisible(lapply(
        X = direction,
        FUN = function(direction) {
            x <- topMarkers(
                data = seurat_all_markers,
                direction = direction
            )
            expect_is(x, "tbl_df")
        }
    ))
})

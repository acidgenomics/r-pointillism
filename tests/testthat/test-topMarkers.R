context("topMarkers")

test_that("grouped_df", {
    x <- topMarkers(all_markers_small)
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

    # Direction
    direction <- formals(topMarkers) %>%
        .[["direction"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(
        X = direction,
        FUN = function(direction) {
            x <- topMarkers(
                data = all_markers_small,
                direction = direction
            )
            expect_is(x, "tbl_df")
        }
    ))
})

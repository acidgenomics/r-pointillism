context("R Markdown")

## FIXME This seems to be hanging on "Importing _setup.R" step at end....

test_that("Seurat per cluster analysis", {
    skeleton <- system.file(
        "rmarkdown",
        "templates",
        "seurat-per-cluster-analysis",
        "skeleton",
        "skeleton.Rmd",
        package = .pkgName,
        mustWork = TRUE
    )
    input <- tempfile(
        pattern = "render",
        tmpdir = tempdir(),
        fileext = ".Rmd"
    )
    file.copy(from = skeleton, to = input, overwrite = TRUE)
    out <- rmarkdown::render(
        input = input,
        output_format = "html_document",
        clean = TRUE,
        params = list(
            "seurat_file" = system.file(
                "data",
                "Seurat.rda",
                package = "AcidTest",
                mustWork = TRUE
            )
        ),
        quiet = TRUE
    )
    expect_true(file.exists(out))
})

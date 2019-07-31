data(
    Seurat,
    package = "acidtest",
    envir = environment()
)
data(
    seuratAllMarkers,
    seuratKnownMarkers,
    package = "pointillism",
    envir = environment()
)

seurat <- Seurat
rm(Seurat)

genes <- head(rownames(seurat))
sce <- as(seurat, "SingleCellExperiment")
objects <- list(
    SingleCellExperiment = sce,
    Seurat = seurat
)

## nolint start
group_vars <- dplyr::group_vars
## nolint end

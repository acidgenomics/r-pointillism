data(
    SingleCellExperiment_Seurat,
    Seurat,
    package = "acidtest",
    envir = environment()
)
sce <- SingleCellExperiment_Seurat
seurat <- Seurat

## FIXME Rename these datasets to camel case.
data(
    seurat_all_markers,
    seurat_known_markers,
    package = "pointillism",
    envir = environment()
)

objects <- list(
    SingleCellExperiment = sce,
    seurat = seurat
)
genes <- head(rownames(seurat))

## nolint start
group_vars <- dplyr::group_vars
## nolint end

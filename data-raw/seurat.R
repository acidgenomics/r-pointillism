# Seurat example data
# 2018-10-21

# # Restrict to 1 MB.
# Use `pryr::object_size()` instead of `utils::object.size()`.
library(pryr)
limit <- structure(1e6, class = "object_size")

# Check and make sure Python umap-learn is accessible to run UMAP.
# We're using this in the `Seurat::RunUMAP()` call below.
# Set `RETICULATE_PYTHON` to conda python binary in `~/.Renviron`.
# This is not working consistently for me on Linux.
library(reticulate)
stopifnot(identical(basename(Sys.getenv("RETICULATE_PYTHON")), "python"))
stopifnot(py_module_available(module = "umap"))

library(splatter)
library(Seurat)
library(Matrix)
library(tidyverse)

library(bcbioSingleCell)
data(pbmc_small, package = "Seurat")
object_size(pbmc_small)
stopifnot(object_size(pbmc_small) < limit)

# seurat_small =================================================================
seurat_small <- pbmc_small %>%
    RunUMAP() %>%
    runZinbwave()
object_size(seurat_small)
stopifnot(object_size(seurat_small) < limit)
validObject(seurat_small)

# `Seurat::pbmc_small` gene symbols map to GRCh37.
gr <- makeGRangesFromEnsembl("Homo sapiens", genomeBuild = "GRCh37")
x <- rownames(seurat_small)
table <- gr$geneName %>%
    as.character() %>%
    make.unique()
names(gr) <- table
stopifnot(all(x %in% table))
which <- match(x = x, table = table)
gr <- gr[which]
rowRanges(seurat_small) <- gr

# all_markers_small ============================================================
markers <- FindAllMarkers(seurat_small)
ranges <- rowRanges(seurat_small)
all_markers_small <- MarkersPerCluster(object = markers, ranges = ranges)

# known_markers_small ==========================================================
gene2symbol <- Gene2Symbol(seurat_small)
data <- DataFrame(
    cellType = as.factor(paste("cell_type", seq_len(2L), sep = "_")),
    geneID = head(gene2symbol[["geneID"]], n = 2L)
)
known_markers_small <- CellTypeMarkers(data, gene2symbol = gene2symbol)
# Write out an example CSV that we can use to test `CellTypeMarkers()`.
export(
    x = do.call(rbind, known_markers_small),
    file = file.path("inst", "extdata", "cell_type_markers.csv")
)

# Save =========================================================================
usethis::use_data(
    seurat_small,
    all_markers_small,
    known_markers_small,
    compress = "xz",
    overwrite = TRUE
)

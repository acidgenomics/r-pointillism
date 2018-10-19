# Seurat example data
# 2018-10-18

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
gr <- makeGRangesFromEnsembl("Homo sapiens", build = "GRCh37")
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
all_markers_small <- SeuratMarkers(
    markers = FindAllMarkers(seurat_small),
    ranges = rowRanges(seurat_small)
)

# known_markers_small ==========================================================
all <- all_markers_small
known <- new(
    Class = "CellTypeMarkers",
    DataFrame(
        cellType = as.factor(paste("cell_type", seq_len(2L), sep = "_")),
        geneID = as.character(head(all$ranges$geneID, n = 2L)),
        geneName = as.character(head(all$ranges$geneName, n = 2L))
    ),
    metadata = list(
        version = packageVersion("pointillism"),
        organism = "Homo sapiens",
        genomeBuild = "GRCh37",
        ensemblRelease = 75L,
        date = Sys.Date()
    )
)
# Write out an example CSV that we can use to test `CellTypeMarkers()`.
export(
    x = known,
    file = file.path("inst", "extdata", "cell_type_markers.csv")
)
known_markers_small <- knownMarkers(all = all, known = known)

# Save =========================================================================
usethis::use_data(
    seurat_small,
    all_markers_small,
    known_markers_small,
    compress = "xz",
    overwrite = TRUE
)

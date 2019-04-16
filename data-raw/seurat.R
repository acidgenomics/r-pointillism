# Seurat example data
# 2019-04-16

library(pryr)
library(reticulate)
library(Matrix)
library(splatter)
library(Seurat)  # 3.0
library(bcbioSingleCell)
library(tidyverse)

# Check and make sure Python umap-learn is accessible to run UMAP.
# We're using this in the `Seurat::RunUMAP()` call below.
#
# umap-learn via reticulate doesn't work well with conda.
# Set up a Python 3 virtual environment instead.
# https://github.com/satijalab/seurat/issues/486
#
# Can also set `RETICULATE_PYTHON` to python binary in `~/.Renviron`.

virtualenv_list()
# [1] "reticulate"

# source ~/.virtualenvs/reticulate/bin/activate
use_virtualenv(virtualenv = "reticulate", required = TRUE)

py_config()

# [Azure VM]
# python:         /home/mike/.virtualenvs/reticulate/bin/python
# libpython:      /usr/local/lib/libpython3.7m.so
# pythonhome:     /usr/local:/usr/local
# version:        3.7.3 (default, Apr 11 2019, 16:32:38)  [GCC 4.8.5 20150623 (Red Hat 4.8.5-36)]
# numpy:          /home/mike/.virtualenvs/reticulate/lib/python3.7/site-packages/numpy
# numpy_version:  1.16.2
#
# python versions found:
#  /usr/bin/python
#  /usr/local/bin/python3
#  /home/mike/.virtualenvs/reticulate/bin/python

# [macOS]
# python:         /usr/local/opt/python/libexec/bin/python
# libpython:      /usr/local/opt/python/Frameworks/Python.framework/Versions/3.7/lib/python3.7/config-3.7m-darwin/libpython3.7.dylib
# pythonhome:     /usr/local/opt/python/Frameworks/Python.framework/Versions/3.7:/usr/local/opt/python/Frameworks/Python.framework/Versions/3.7
# version:        3.7.2 (default, Feb 12 2019, 08:15:36)  [Clang 10.0.0 (clang-1000.11.45.5)]
# numpy:          /usr/local/lib/python3.7/site-packages/numpy
# numpy_version:  1.16.2
#
# python versions found:
#  /usr/local/opt/python/libexec/bin/python
#  /usr/bin/python
#  /usr/local/bin/python
#  /usr/local/bin/python3
#  /Users/mike/anaconda3/bin/python
#  /Users/mike/.virtualenvs/reticulate/bin/python

stopifnot(py_module_available(module = "umap"))

# # Restrict object size to 1 MB.
# Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(1e6, class = "object_size")

# seurat =======================================================================
data(pbmc_small, package = "Seurat")
object_size(pbmc_small)
stopifnot(object_size(pbmc_small) < limit)

# The Seurat wiki describes the changes in v3.0.
# https://github.com/satijalab/seurat/wiki
seurat <- pbmc_small

# Add UMAP dimensional reduction to example object.
# Alternatively, can use `features` here instead.
seurat <- RunUMAP(seurat, dims = seq_len(10L))

object_size(seurat)
stopifnot(object_size(seurat) < limit)
validObject(seurat)

# The pbmc_small gene symbols map to GRCh37.
gr <- makeGRangesFromEnsembl(organism = "Homo sapiens", genomeBuild = "GRCh37")
x <- rownames(seurat)
table <- make.unique(as.character(gr$geneName))
names(gr) <- table
stopifnot(all(x %in% table))
which <- match(x = x, table = table)
gr <- gr[which]
# Note that `rowRanges` method is defined by pointillism.
rowRanges(seurat) <- gr

# seurat_all_markers ===========================================================
seurat_all_markers <- SeuratMarkersPerCluster(
    object = FindAllMarkers(seurat),
    ranges = rowRanges(seurat)
)

# known_markers_small ==========================================================
data(cell_type_markers)
seurat_known_markers <- KnownMarkers(
    markers = seurat_all_markers,
    known = cell_type_markers$homoSapiens
)
export(
    x = seurat_known_markers,
    file = file.path("inst", "extdata", "cell_type_markers.csv")
)

# Save =========================================================================
usethis::use_data(
    seurat,
    seurat_all_markers,
    seurat_known_markers,
    compress = "xz",
    overwrite = TRUE
)

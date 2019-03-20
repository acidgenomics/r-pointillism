# Seurat example data
# 2019-03-20

library(reticulate)
library(Seurat)

library(pryr)
library(splatter)
library(Matrix)
library(tidyverse)
library(bcbioSingleCell)

# Check and make sure Python umap-learn is accessible to run UMAP.
# We're using this in the `Seurat::RunUMAP()` call below.
# Set `RETICULATE_PYTHON` to pip virtualenv python binary in `~/.Renviron`.

# umap-learn via reticulate doesn't work well with conda.
# https://github.com/satijalab/seurat/issues/486

# ~/.virtualenvs/reticulate/bin/activate
# use_virtualenv(virtualenv = "reticulate")

use_virtualenv(virtualenv = "reticulate")

py_config()
# python:         /home/mike/.virtualenvs/reticulate/bin/python
# libpython:      /usr/lib64/python2.7/config/libpython2.7.so
# pythonhome:     /usr:/usr
# virtualenv:     /home/mike/.virtualenvs/reticulate/bin/activate_this.py
# version:        2.7.5 (default, Sep 12 2018, 05:31:16)  [GCC 4.8.5 20150623 (Red Hat 4.8.5-36)]
# numpy:          /home/mike/.virtualenvs/reticulate/lib/python2.7/site-packages/numpy
# numpy_version:  1.16.2

stopifnot(py_module_available(module = "umap"))

# # Restrict object size to 1 MB.
# Use `pryr::object_size` instead of `utils::object.size`.
limit <- structure(1e6, class = "object_size")

data(pbmc_small, package = "Seurat")
object_size(pbmc_small)
stopifnot(object_size(pbmc_small) < limit)

# seurat_small =================================================================
seurat_small <- pbmc_small
# Note that this step requires umap-learn via reticulate.

# Python 2.7 issue?
# Error in py_call_impl(callable, dots$args, dots$keywords) :
#  Evaluation error: Required version of NumPy not available: installation of Numpy >= 1.6 not found.
# Calls: RunUMAP -> <Anonymous> -> py_call_impl

# Error in py_call_impl(callable, dots$args, dots$keywords) :
#   TypingError: Failed in nopython mode pipeline (step: nopython frontend)
# Unknown attribute 'shape' of type none
#
# File "../../../../../home/mike/.virtualenvs/reticulate/lib/python2.7/site-packages/umap/umap_.py", line 88:
# def smooth_knn_dist(distances, k, n_iter=64, local_connectivity=1.0, bandwidth=1.0):
#     <source elided>
#     target = np.log2(k) * bandwidth
#     rho = np.zeros(distances.shape[0])
#     ^
#
# [1] During: typing of get attribute at /home/mike/.virtualenvs/reticulate/lib/python2.7/site-packages/umap/umap_.py (88)
#
# File "../../../../../home/mike/.virtualenvs/reticulate/lib/python2.7/site-packages/umap/umap_.py", line 88:
# def smooth_knn_dist(distances, k, n_iter=64, local_connectivity=1.0, bandwidth=1.0):
#     <source elided>
#     target = np.log2(k) * bandwidth
#     rho = np.zeros(distances.shape[0])
#     ^
#
# This is not usually a problem with Numba itself but instead often caused by
# the use of unsupported features or an issue in resolving types.
#
# To see Python/NumPy features supported by the latest release of Numba visit:
# http://numba.pydata.org/numba-doc/dev/reference/pysupported.html
# and
# http://numba.pydata.org/numba-doc/dev/reference/numpysupported.html
#
# For more information about typing errors and how to debug them visit:
# http://numba.pydata.org/numba-doc/latest/user/troubleshoot.html#my-code-doesn-t-compile
#
# If you think your code should work with Numba, please report the error message
# and traceback, along with a minimal reproducer at:
# https://github.com/numba/numba/issues/new
#
# Detailed traceback:
#   File "/home/mike/.virtualenvs/reticulate/lib/python2.7/site-packages/umap/umap_.py", line 1566, in fit_transform
#     self.fit(X, y)
#   File "/home/mike/.virtualenvs/reticulate/lib/python2.7/site-packages/umap/umap_.py", line 1398, in fit
#     self.verbose,
#   File "/home/mike/.virtualenvs/reticulate/lib/python2.7/site-packages/numba/dispatcher.py", line 350, in _compile_for_args
#     error_rewrite(e, 'typing')
#   File "/home/mike/.virtualenvs/reticulate/lib/python2.7/site-packages/numba/dispatcher.py", line 317, in error_rewrite
#     reraise(type(e), e, None)
#   File "<string>", line 2, in reraise
# Calls: RunUMAP -> <Anonymous> -> py_call_impl

seurat_small <- RunUMAP(seurat_small)
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
all_markers_small <- SeuratMarkersPerCluster(object = markers, ranges = ranges)

# known_markers_small ==========================================================
known_markers_small <- KnownMarkers(
    markers = all_markers_small,
    known = cell_type_markers$homoSapiens
)
export(
    x = known_markers_small,
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

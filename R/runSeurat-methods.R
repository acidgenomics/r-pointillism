#' Run Seurat
#'
#' @name runSeurat
#' @note Updated 2020-06-26.
#'
#' @section reticulate:
#'
#' A Python virtual environment (e.g. "r-reticulate") with umap-learn installed
#' is required when `umapMethod` is set to "umap-learn" (default).
#'
#' The reticulate package supports conda or virtual environments managed by
#' venv/virtualenv. I've had better luck with using using virtual environments
#' rather than attempting to use conda with this package.
#'
#' Check `Sys.getenv("WORKON_HOME")` path, which is the current approach used
#' in `Rprofile.site`.
#'
#' @inheritParams AcidRoxygen::params
#' @param regressCellCycle `character(1)`.
#'   - `"s-g2m-diff"`: Calculate the difference between S and G2/M phases and
#'     use that to regress. See `CC.Difference` metric in Seurat vignette.
#'   - `"yes"`: Regress out any effects of both S and G2/M phase variable.
#'     Refer to `"S.Score"` and `"G2M.Score"` metrics in Seurat vignette.
#'   - `"no"`: Don't calculate cell-cycle scoring and don't regress.
#'
#'   Refer to the Seurat cell-cycle regression vignette for details.
#' @param varsToRegress `character` or `NULL`.
#'   Unwanted sources of variance to regress. Note that when `regressCellCycle`
#'   is not `"no"`, then the corresponding cell-cycle variables are added
#'   automatically. Passes to [Seurat::ScaleData] internally.
#' @param dims `"auto"` or `integer`.
#'   Dimensions of reduction to use as input for shared nearest neighbor (SNN)
#'   graph construction. When set to "auto" (default), the elbow point is
#'   calculated internally. See [plotPCElbow()] for details. Passes to
#'   [Seurat::FindNeighbors()] and [Seurat::RunUMAP()] internally.
#' @param resolution `numeric`.
#'   Resolutions to calculate for clustering.
#'   Passes to [Seurat::FindClusters()] internally.
#' @param tsneMethod `character(1)`.
#'   t-SNE algorithm passed to [Seurat::RunTSNE()].
#'   Using the current default in Seurat.
#' @param umapMethod `character(1)`.
#'   UMAP algorithm passed to [Seurat::RunUMAP()].
#'
#'   Currently supported:
#'   - `"uwot"`, changed to default in Seurat 3.
#'     Note that this sets `metric = "cosine"` automatically.
#'   - `"umap-learn"`, which requires reticulate.
#'     Note that this sets `metric = "correlation"` automatically.
#' @param virtualenv `character(1)`.
#'   Python virtual environment name.
#'   Only evaluated when `umapMethod = "umap-learn"`.
#'   See reticulate section for details.
#' @param workers `"auto"`, `integer(1)`, or `NULL`.
#'   Disable parallelization with future by setting to `NULL`.
#'
#' @return `Seurat`.
#'
#' @seealso
#' - https://github.com/satijalab/seurat/wiki
#' - https://satijalab.org/seurat/essential_commands.html
#' - https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
#'   https://satijalab.org/seurat/v3.0/future_vignette.html
#' - https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
NULL



## Updated 2020-06-26.
`runSeurat,Seurat` <-  # nolint
    function(
        object,
        regressCellCycle = c("s-g2m-diff", "yes", "no"),
        varsToRegress = c("nCount_RNA", "mitoRatio"),
        dims = "auto",
        resolution = seq(from = 0.2, to = 1.2, by = 0.2),
        tsneMethod = "Rtsne",
        umapMethod = "uwot",
        virtualenv = "r-reticulate",
        workers = "auto"
    ) {
        assert(
            requireNamespace("future", quietly = TRUE),
            requireNamespace("Seurat", quietly = TRUE),
            isCharacter(varsToRegress, nullOK = TRUE),
            identical(dims, "auto") || is.numeric(dims),
            is.numeric(resolution),
            isString(tsneMethod),
            isString(umapMethod),
            isString(virtualenv),
            identical(workers, "auto") || isInt(workers, nullOK = TRUE)
        )
        regressCellCycle <- match.arg(regressCellCycle)
        ## Parallelization (via future) ----------------------------------------
        ## Note that Seurat currently uses future package for parallelization.
        ## Multiprocess is currently unstable in RStudio and disabled.
        if (
            isTRUE(future::supportsMulticore()) &&
            !is.null(workers)
        ) {
            if (identical(workers, "auto")) {
                workers <- max(getOption(x = "mc.cores", default = 1L), 1L)
            }
            assert(isInt(workers))
            alert(sprintf(
                "Enabling {.pkg %s} multiprocess with %d workers.",
                "future", workers
            ))
            future::plan("multiprocess", workers = workers)
        }
        ## Pre-processing ------------------------------------------------------
        alert(sprintf(
            "{.pkg %s}::{.fun %s}",
            "Seurat", "NormalizeData"
        ))
        object <- Seurat::NormalizeData(object)
        alert(sprintf(
            "{.pkg %s}::{.fun %s}",
            "Seurat", "FindVariableFeatures"
        ))
        object <- Seurat::FindVariableFeatures(object)
        ## Cell-cycle regression -----------------------------------------------
        ## https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
        if (!identical(regressCellCycle, "no")) {
            env <- new.env()
            data(
                "cell_cycle_markers_list",
                package = "pointillism",
                envir = env
            )
            ccm <- env[["cell_cycle_markers_list"]]
            organism <- camelCase(organism(object), strict = TRUE)
            if (!isSubset(organism, names(ccm))) {
                abort(sprintf(
                    fmt = paste(
                        "Failed to obtain cell-cycle markers.",
                        "{.val %s} is not currently supported.",
                        "Please file an issue on GitHub: {.url %s}.",
                        sep = "\n"
                    ),
                    organism,
                    "https://github.com/acidgenomics/pointillism/issues"
                ))
            }
            ccm <- ccm[[organism]]
            assert(is(ccm, "CellCycleMarkers"))
            g2mGenes <- as.character(ccm[["g2m"]][["geneName"]])
            sGenes <- as.character(ccm[["s"]][["geneName"]])
            alert(sprintf(
                "{.pkg %s}::{.fun %s}",
                "Seurat", "CellCycleScoring"
            ))
            object <- Seurat::CellCycleScoring(
                object = object,
                s.features = sGenes,
                g2m.features = g2mGenes
            )
        }
        ## Update the variables to regress, including cell cycle.
        if (identical(regressCellCycle, "s-g2m-diff")) {
            ## Note that Seurat uses non-standard '$' and '[[' methods.
            object$CC.Difference <- object$S.Score - object$G2M.Score  # nolint
            varsToRegress <- c("CC.Difference", varsToRegress)
        } else if (identical(regressCellCycle, "yes")) {
            varsToRegress <- c("S.Score", "G2M.Score", varsToRegress)
        }
        ## Scaling and projection ----------------------------------------------
        alert(sprintf(
            "{.pkg %s}::{.fun %s}",
            "Seurat", "ScaleData"
        ))
        dl(c("varsToRegress" = toInlineString(varsToRegress)))
        ## Scaling all features is very slow for large datasets.
        ## Current default in Seurat scales variable features only.
        ## All features:
        ## > features = rownames(object)
        ## Variable features only:
        ## > features = Seurat::VariableFeatures(object)
        object <- Seurat::ScaleData(
            object = object,
            features = rownames(object),
            vars.to.regress = varsToRegress
        )
        alert(sprintf(
            "{.pkg %s}::{.fun %s}",
            "Seurat", "RunPCA"
        ))
        object <- Seurat::RunPCA(object)
        if (identical(dims, "auto")) {
            alert(sprintf("{.fun %s}", "plotElbow"))
            p <- plotPCElbow(object)
            elbow <- attr(p, "elbow")
            dims <- seq(from = 1L, to = elbow, by = 1L)
        }
        ## Clustering ----------------------------------------------------------
        alert(sprintf(
            "{.pkg %s}::{.fun %s}",
            "Seurat", "FindNeighbors"
        ))
        alertInfo(sprintf("Using %d dims,", length(dims)))
        object <- Seurat::FindNeighbors(object, dims = dims)
        alert(sprintf(
            "{.pkg %s}::{.fun %s}",
            "Seurat", "FindClusters"
        ))
        dl(c("resolution" = as.character(resolution)))
        object <- Seurat::FindClusters(object, resolution = resolution)
        ## tSNE / UMAP ---------------------------------------------------------
        alert(sprintf(
            "{.pkg %s}::{.fun %s}",
            "Seurat", "RunTSNE"
        ))
        dl(c("method" = tsneMethod))
        object <- Seurat::RunTSNE(
            object = object,
            tsne.method = tsneMethod
        )
        if (identical(umapMethod, "umap-learn")) {
            assert(requireNamespace("reticulate", quietly = TRUE))
            reticulate::use_virtualenv(virtualenv = virtualenv, required = TRUE)
            assert(reticulate::py_module_available(module = "umap"))
        }
        alert(sprintf(
            "{.pkg %s}::{.fun %s}",
            "Seurat", "RunUMAP"
        ))
        metric <- switch(
            EXPR = umapMethod,
            "umap-learn" = "correlation",
            "uwot" = "cosine"
        )
        dl(c(
            "method" = umapMethod,
            "metric" = metric
        ))
        alertInfo(paste("Using", length(dims), "dims."))
        ## The default method for RunUMAP has changed from calling Python UMAP
        ## via reticulate to the R-native UWOT using the cosine metric. To use
        ## Python UMAP via reticulate, set umap.method to 'umap-learn' and
        ## metric to 'correlation'.
        object <- Seurat::RunUMAP(
            object = object,
            umap.method = umapMethod,
            metric = metric,
            dims = dims
        )
        alertSuccess("Seurat run was successful.")
        object
    }



## Updated 2020-06-26.
`runSeurat,SCE` <-  # nolint
    function(object, ...) {
        object <- as(object, "SingleCellExperiment")
        object <- convertGenesToSymbols(object)
        object <- as(object, "Seurat")
        runSeurat(object, ...)
    }



#' @rdname runSeurat
#' @export
setMethod(
    f = "runSeurat",
    signature = signature(object = "Seurat"),
    definition = `runSeurat,Seurat`
)

#' @rdname runSeurat
#' @export
setMethod(
    f = "runSeurat",
    signature = signature(object = "SingleCellExperiment"),
    definition = `runSeurat,SCE`
)

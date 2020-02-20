#' Run Seurat
#'
#' @export
#' @note Updated 2020-02-20.
#'
#' @section reticulate:
#'
#' Check `Sys.getenv("WORKON_HOME")` path, which is the current approach used
#' in `Rprofile.site`.
#'
#' @seealso
#' - https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#' - https://github.com/satijalab/seurat/wiki
#' - https://satijalab.org/seurat/essential_commands.html
#' - https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
#' - Parallelization with future
#'   https://satijalab.org/seurat/v3.0/future_vignette.html
runSeurat <- function(
    object,
    regressCellCycle = c("s-g2m-diff", "yes", "no"),
    varsToRegress = c("nCount_RNA", "mitoRatio"),
    dims = "auto",
    resolution = seq(from = 0.2, to = 1.2, by = 0.2),
    tsneMethod = "Rtsne",
    umapMethod = "umap-learn",
    virtualenv = "r-reticulate",
    workers = "auto"
) {
    assert(
        requireNamespace("future", quietly = TRUE),
        requireNamespace("Seurat", quietly = TRUE),
        is(object, "Seurat"),
        isCharacter(varsToRegress, nullOK = TRUE),
        identical(dims, "auto") || is.numeric(dims),
        is.numeric(resolution),
        isString(tsneMethod),
        isString(umapMethod),
        isString(virtualenv),
        identical(dims, "auto") || isInt(workesr)
    )
    regressCellCycle <- match.arg(regressCellCycle)

    ## Parallelization (via future) --------------------------------------------
    ## Note that Seurat currently uses future package for parallelization.
    ## Multiprocess is currently unstable in RStudio and disabled.
    if (isTRUE(future::supportsMulticore())) {
        if (identical(workers, "auto")) {
            workers <- max(getOption("mc.cores"), 1L)
        }
        assert(isInt(workers))
        cli_alert(paste(
            "Enabling {.pkg future} multiprocess with",
            workers, "workers."
        ))
        future::plan("multiprocess", workers = workers)
    }

    ## Pre-processing ----------------------------------------------------------
    cli_alert("{.pkg Seurat}::{.fun NormalizeData}")
    object <- Seurat::NormalizeData(object)

    cli_alert("{.pkg Seurat}::{.fun FindVariableFeatures}")
    object <- Seurat::FindVariableFeatures(object)

    ## Cell-cycle regression ---------------------------------------------------
    ## https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html

    if (!identical(regressCellCycle, "no")) {
        ## FIXME This step will error if we subset the Seurat object...
        organism <- organism(object)
        data(
            cell_cycle_markers_list,
            package = "pointillism",
            envir = environment()
        )
        ccm <- cell_cycle_markers_list[[camelCase(organism)]]
        assert(is(ccm, "CellCycleMarkers"))
        g2mGenes <- as.character(ccm[["g2m"]][["geneName"]])
        sGenes <- as.character(ccm[["s"]][["geneName"]])
        cli_alert("{.pkg Seurat}::{.fun CellCycleScoring}")
        object <- Seurat::CellCycleScoring(
            object = object,
            s.features = sGenes,
            g2m.features = g2mGenes
        )
    }

    ## Update the variables to regress, including cell cycle.
    if (identical(regressCellCycle, "s-g2m-diff")) {
        object[["CC.Difference"]] <- object$S.Score - object$G2M.Score
        varsToRegress <- c("CC.Difference", varsToRegress)
    } else if (identical(regressCellCycle, "yes")) {
        varsToRegress <- c("S.Score", "G2M.Score", varsToRegress)
    }

    ## Scaling and projection --------------------------------------------------
    cli_alert("{.pkg Seurat}::{.fun ScaleData}")
    cli_dl(c(varsToRegress = toString(varsToRegress)))
    ## This step is CPU intensive for large datasets.
    object <- Seurat::ScaleData(
        object = object,
        vars.to.regress = varsToRegress,
        features = rownames(object)
    )

    cli_alert("{.pkg Seurat}::{.fun RunPCA}")
    object <- Seurat::RunPCA(
        object = object,
        features = Seurat::VariableFeatures(object)
    )

    if (identical(dims, "auto")) {
        cli_alert("{.fun plotElbow}")
        dims <- plotElbow(object)
    }

    ## Clustering --------------------------------------------------------------
    cli_alert("{.pkg Seurat}::{.fun FindNeighbors}")
    object <- Seurat::FindNeighbors(object, dims = dims)

    cli_alert("{.pkg Seurat}::{.fun FindClusters}")
    object <- Seurat::FindClusters(object, resolution = resolution)

    if (identical(umapMethod, "umap-learn")) {
        assert(requireNamespace("reticulate", quietly = TRUE))
        reticulate::use_virtualenv(virtualenv = virtualenv, required = TRUE)
        assert(reticulate::py_module_available(module = "umap"))
    }

    ## tSNE / UMAP -------------------------------------------------------------
    cli_alert("{.pkg Seurat}::{.fun RunTSNE}")
    object <- Seurat::RunTSNE(
        object = object,
        tsne.method = tsneMethod
    )

    cli_alert("{.pkg Seurat}::{.fun RunUMAP}")
    object <- Seurat::RunUMAP(
        object = object,
        umap.method = umapMethod,
        dims = dims
    )

    object
}

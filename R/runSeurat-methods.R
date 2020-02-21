#' Run Seurat
#'
#' @name runSeurat
#' @note Updated 2020-02-21.
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
NULL



## Updated 2020-02-21.
`runSeurat,Seurat` <-  # nolint
    function(
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
            isCharacter(varsToRegress, nullOK = TRUE),
            identical(dims, "auto") || is.numeric(dims),
            is.numeric(resolution),
            isString(tsneMethod),
            isString(umapMethod),
            isString(virtualenv),
            identical(workers, "auto") || isInt(workers)
        )
        regressCellCycle <- match.arg(regressCellCycle)

        ## Parallelization (via future) ----------------------------------------
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

        ## Pre-processing ------------------------------------------------------
        cli_alert("{.pkg Seurat}::{.fun NormalizeData}")
        object <- Seurat::NormalizeData(object)

        cli_alert("{.pkg Seurat}::{.fun FindVariableFeatures}")
        object <- Seurat::FindVariableFeatures(object)

        ## Cell-cycle regression -----------------------------------------------
        ## https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html

        if (!identical(regressCellCycle, "no")) {
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
            ## Note that Seurat uses non-standard '$' and '[[' methods.
            object$CC.Difference <- object$S.Score - object$G2M.Score  # nolint
            varsToRegress <- c("CC.Difference", varsToRegress)
        } else if (identical(regressCellCycle, "yes")) {
            varsToRegress <- c("S.Score", "G2M.Score", varsToRegress)
        }

        ## Scaling and projection ----------------------------------------------
        cli_alert("{.pkg Seurat}::{.fun ScaleData}")
        cli_dl(c(varsToRegress = toString(varsToRegress)))
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

        cli_alert("{.pkg Seurat}::{.fun RunPCA}")
        object <- Seurat::RunPCA(object)

        if (identical(dims, "auto")) {
            cli_alert("{.fun plotElbow}")
            p <- plotPCElbow(object)
            elbow <- attr(p, "elbow")
            dims <- seq(from = 1L, to = elbow, by = 1L)
        }

        ## Clustering ----------------------------------------------------------
        cli_alert("{.pkg Seurat}::{.fun FindNeighbors}")
        cli_alert_info(paste("Using", length(dims), "dims."))
        object <- Seurat::FindNeighbors(object, dims = dims)

        cli_alert("{.pkg Seurat}::{.fun FindClusters}")
        cli_dl(c(resolution = deparse(resolution)))
        object <- Seurat::FindClusters(object, resolution = resolution)

        ## tSNE / UMAP ---------------------------------------------------------
        cli_alert("{.pkg Seurat}::{.fun RunTSNE}")
        cli_dl(c(method = tsneMethod))
        object <- Seurat::RunTSNE(
            object = object,
            tsne.method = tsneMethod
        )

        if (identical(umapMethod, "umap-learn")) {
            assert(requireNamespace("reticulate", quietly = TRUE))
            reticulate::use_virtualenv(virtualenv = virtualenv, required = TRUE)
            assert(reticulate::py_module_available(module = "umap"))
        }

        cli_alert("{.pkg Seurat}::{.fun RunUMAP}")
        cli_dl(c(method = umapMethod))
        cli_alert_info(paste("Using", length(dims), "dims."))
        object <- Seurat::RunUMAP(
            object = object,
            umap.method = umapMethod,
            dims = dims
        )

        cli_alert_success("Seurat run was successful.")
        object
    }



#' @rdname runSeurat
#' @export
setMethod(
    f = "runSeurat",
    signature = signature("Seurat"),
    definition = `runSeurat,Seurat`
)



## Updated 2020-02-20.
`runSeurat,SingleCellExperiment` <-  # nolint
    function(object, ...) {
        runSeurat(
            object = as(object, "Seurat"),
            ...
        )
    }



#' @rdname runSeurat
#' @export
setMethod(
    f = "runSeurat",
    signature = signature("SingleCellExperiment"),
    definition = `runSeurat,SingleCellExperiment`
)

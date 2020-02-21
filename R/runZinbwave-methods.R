#' Run ZINB-WaVE
#'
#' @name runZinbwave
#' @inherit zinbwave::zinbwave title description params
#' @note Updated 2020-02-21.
#'
#' @details
#' [zinbwave][] will calculate `normalizedValues` and `weights` matrices, which
#' will be defined in the `assays` slot. If these values have already been
#' calculated, the object will be returned unmodified, unless `recalculate =
#' TRUE`.
#'
#' [zinbwave]: https://bioconductor.org/packages/release/bioc/html/zinbwave.html
#'
#' @section Parallelization:
#'
#' `zinbwave()` defaults to using `BiocParallel::bpparam` to register the number
#' of cores to use (`bpparam` argument), but we've found that this works
#' inconsistently across installations. Currently, we recommend using
#' `BiocParallel::SerialParam` by default, which will run in serial, using a
#' single core. On Linux or macOS, `BiocParallel::MulticoreParam` should work
#' to run with multiple cores.
#'
#' @inheritParams acidroxygen::params
#' @inheritParams zinbwave::zinbFit
#' @inheritParams zinbwave::zinbModel
#'
#' @param recalculate `logical`. Force recalculation of weights.
#' @param verbose `logical`. Run zinbwave in verbose mode (for debugging).
#'
#' @return Modified object.
#'   Weights added to `assays`.
#'
#' @seealso `zinbwave::zinbwave`.
#'
#' @examples
#' data(SingleCellExperiment, package = "acidtest")
#' Y <- SingleCellExperiment
#' Y <- nonzeroRowsAndCols(Y)
#' Y <- runZinbwave(Y)
#' class(weights(Y))
NULL



## Updated 2020-02-21.
`runZinbwave,SingleCellExperiment` <-  # nolint
    function(
        Y,  # nolint
        K = 0L,  # nolint
        epsilon = 1e12,
        # Use serial by default, for cross-platform compatibility.
        BPPARAM,  # nolint
        recalculate = FALSE,
        verbose = FALSE
    ) {
        assert(
            requireNamespace("zinbwave", quietly = TRUE),
            is(Y, "SingleCellExperiment"),
            isNumber(K),
            isNumber(epsilon),
            .isBPPARAM(BPPARAM),
            isFlag(recalculate),
            isFlag(verbose)
        )
        cli_alert("Running zinbwave.")
        # Early return if weights are calculated.
        weights <- tryCatch(
            expr = weights(Y),
            error = function(e) NULL
        )
        if (is.matrix(weights) && !isTRUE(recalculate)) {
            message("Object already contains pre-calculated weights.")
            return(Y)
        }

        # BiocParallel ---------------------------------------------------------
        # Use a progress bar (only applies to multicore).
        bpprogressbar(BPPARAM) <- TRUE
        # Inform the user whether running in parallel or serial.
        bpparamInfo <- capture.output(BPPARAM)
        message(paste(
            "BiocParallel param registered.",
            paste(bpparamInfo, collapse = "\n"),
            sep = "\n"
        ))

        # Prepare Y object for zinbwave ----------------------------------------
        # We're returnining original object class, with modified assays.
        # Note that `zinbwave` otherwise returns `SingleCellExperiment`.
        object <- Y
        # Coerce to SingleCellExperiment, for consistency.
        Y <- as(Y, "SingleCellExperiment")  # nolint
        # zinbFit doesn't currently support sparse counts.
        # Ensure they are coerced to a dense matrix.
        # Keep an original copy in case they're sparse, and reslot.
        counts(Y) <- as.matrix(counts(Y))

        # Fit a ZINB regression model ------------------------------------------
        message("zinbFit(): Fitting a ZINB regression model.")
        message(paste(
            "CPU time used:",
            printString(system.time({
                # Wrapping here to disable progress bars.
                invisible(capture.output(
                    fittedModel <- zinbwave::zinbFit(
                        Y = Y,
                        K = K,
                        epsilon = epsilon,
                        BPPARAM = BPPARAM,
                        verbose = verbose
                    )
                ))
            })),
            sep = "\n"
        ))
        message(paste(
            capture.output(print(fittedModel)),
            collapse = "\n"
        ))

        # zinbwave -------------------------------------------------------------
        message(paste(
            "zinbwave(): Performing dimensionality reduction using a",
            "ZINB regression model with gene and cell-level covariates."
        ))
        message(paste(
            "CPU time used:",
            printString(system.time({
                # Wrapping here to disable progress bars.
                invisible(capture.output(
                    zinb <- zinbwave::zinbwave(
                        Y = Y,
                        fitted_model = fittedModel,
                        K = K,
                        epsilon = epsilon,
                        BPPARAM = BPPARAM,
                        verbose = verbose
                    )
                ))
            })),
            sep = "\n"
        ))

        # Return ---------------------------------------------------------------
        # Re-slot original raw counts, in case they are sparse.
        counts <- counts(object)
        assert(identical(dimnames(zinb), dimnames(counts)))
        assays(object) <- assays(zinb)
        counts(object) <- counts
        metadata(object)[["weights"]] <- "zinbwave"
        object
    }

args <- "BPPARAM"
formals(`runZinbwave,SingleCellExperiment`)[args] <- .formalsList[args]



#' @rdname runZinbwave
#' @export
setMethod(
    f = "runZinbwave",
    signature = signature(Y = "SingleCellExperiment"),
    definition = `runZinbwave,SingleCellExperiment`
)



## > ## Updated 2020-02-21.
## > `runZinbwave,Seurat` <-  # nolint
## >     function(Y, ...) {  # nolint
## >         zinb <- runZinbwave(Y = as(Y, "SingleCellExperiment"), ...)
## >         weights(Y) <- weights(zinb)
## >         metadata(Y)[["weights"]] <- "zinbwave"
## >         Y
## >     }
## >
## >
## >
## > #' @describeIn runZinbwave Coerces to `SingleCellExperiment` and stashes weights
## > #'   in the `weights()` slot.
## > #' @export
## > setMethod(
## >     f = "runZinbwave",
## >     signature = signature("Seurat"),
## >     definition = `runZinbwave,Seurat`
## > )

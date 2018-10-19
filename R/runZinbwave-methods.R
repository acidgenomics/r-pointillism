# FIXME Add a show progress option.



#' Calculate ZINB-WaVE Weights
#'
#' zinbwave will calculate `normalizedValues` and `weights` matrices, which will
#' be defined in the [assays()] slot. If these values have already been
#' calculated, the object will be returned unmodified, unless
#' `recalculate = TRUE`.
#'
#' @note zinbwave defaults to using [BiocParallel::bpparam()] to register the
#'   number of cores to use (`bpparam` argument), but we've found that this
#'   works inconsistently across installations. Currently, we recommend using
#'   [BiocParallel::SerialParam()] by default, which will run in serial, using a
#'   single core. On Linux or macOS, [BiocParallel::MulticoreParam()] should
#'   work to run with multiple cores.
#'
#' @name runZinbwave
#' @family Differential Expression Functions
#' @include globals.R
#'
#' @inheritParams general
#' @inheritParams zinbwave::zinbwave
#' @inheritParams zinbwave::zinbFit
#' @inheritParams zinbwave::zinbModel
#' @param recalculate `logical`. Force recalculation of weights.
#' @param verbose `logical`. Run zinbwave in verbose model (for debugging).
#' @param ... Passthrough arguments to [zinbwave::zinbwave()].
#'
#' @return Modified object (S4 class must extend `SingleCellExperiment`), with
#'   weights added to [assays()].
#'
#' @seealso [zinbwave::zinbwave()].
#'
#' @examples
#' data(seurat_small)
#' object <- seurat_small
#' ## Example using 100 genes, 100 cells.
#' Y <- as(object, "SingleCellExperiment")
#' Y <- Y[seq_len(100L), seq_len(100L)]
#' zinb <- suppressMessages(runZinbwave(
#'     Y = Y,
#'     BPPARAM = BiocParallel::SerialParam(),
#'     recalculate = TRUE
#' ))
#' print(zinb)
#' class(weights(zinb))
NULL



# SingleCellExperiment =========================================================
.runZinbwave.SCE <-  # nolint
    function(
        Y,  # nolint
        K = 0L,  # nolint
        epsilon = 1e12,
        # Use serial by default, for cross-platform compatibility.
        BPPARAM,  # nolint
        recalculate = FALSE,
        verbose = FALSE,
        ...
    ) {
        message("Running zinbwave.")
        stopifnot(is(Y, "SingleCellExperiment"))

        # Early return if weights are calculated -------------------------------
        weights <- tryCatch(
            expr = weights(Y),
            error = function(e) NULL
        )
        if (is.matrix(weights) && !isTRUE(recalculate)) {
            warning(
                "Object already contains pre-calculated weights.",
                call. = FALSE
            )
            return(Y)
        }

        # Assert checks --------------------------------------------------------
        assert_is_a_number(K)
        assert_is_a_number(epsilon)
        .assertIsBPPARAM(BPPARAM)
        assert_is_a_bool(recalculate)
        assert_is_a_bool(verbose)

        # BiocParallel ---------------------------------------------------------
        # Use a progress bar (only applies to multicore).
        bpprogressbar(BPPARAM) <- TRUE

        # Inform the user whether running in parallel or serial.
        bpparamInfo <- capture.output(BPPARAM)
        message(paste(
            "BiocParallel param registered:",
            paste(bpparamInfo, collapse = "\n"),
            sep = "\n"
        ))

        # Prepare Y object for zinbwave ----------------------------------------
        # We're returnining original object class, with modified assays.
        # Note that `zinbwave()` otherwise returns `SingleCellExperiment`.
        object <- Y

        # Coerce to SingleCellExperiment, for consistency.
        Y <- as(Y, "SingleCellExperiment")  # nolint

        # zinbFit doesn't currently support sparse counts.
        # Ensure they are coerced to a dense matrix.
        # Keep an original copy in case they're sparse, and reslot.
        counts(Y) <- as.matrix(counts(Y))

        # Fit a ZINB regression model ------------------------------------------
        message("Fitting a ZINB regression model.")
        message(paste(
            "CPU time used:",
            printString(system.time({
                fittedModel <- zinbFit(
                    Y = Y,
                    K = K,
                    epsilon = epsilon,
                    BPPARAM = BPPARAM,
                    verbose = verbose
                )
            })),
            sep = "\n"
        ))
        message(paste(
            capture.output(print(fittedModel)),
            collapse = "\n"
        ))

        # zinbwave -------------------------------------------------------------
        message("Running zinbwave.")
        message(paste(
            "CPU time used:",
            printString(system.time({
                zinb <- zinbwave(
                    Y = Y,
                    fitted_model = fittedModel,
                    K = K,
                    epsilon = epsilon,
                    BPPARAM = BPPARAM,
                    verbose = verbose,
                    ...
                )
            })),
            sep = "\n"
        ))

        # Return ---------------------------------------------------------------
        # Re-slot original raw counts, in case they are sparse.
        counts <- counts(object)
        assert_are_identical(dimnames(zinb), dimnames(counts))
        assays(object) <- assays(zinb)
        counts(object) <- counts
        metadata(object)[["weights"]] <- "zinbwave"
        object
    }
formals(.runZinbwave.SCE)[["BPPARAM"]] <- bpparam



#' @rdname runZinbwave
#' @export
setMethod(
    f = "runZinbwave",
    signature = signature(Y = "SingleCellExperiment"),
    definition = .runZinbwave.SCE
)



# seurat =======================================================================
.runZinbwave.seurat <- function(Y, ...) {
    zinb <- runZinbwave(Y = as(Y, "SingleCellExperiment"), ...)
    weights(Y) <- weights(zinb)
    metadata(Y)[["weights"]] <- "zinbwave"
    Y
}



#' @describeIn runZinbwave Coerces to `SingleCellExperiment` and stashes weights
#'   in the [weights()] slot.
#' @export
setMethod(
    f = "runZinbwave",
    signature = signature("seurat"),
    definition = .runZinbwave.seurat
)

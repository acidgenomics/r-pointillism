## Require zero weights to be calculated, showing a more informative error.
.weights <- function(object) {
    msg <- paste(
        "Object does not contain zero weights in `weights` slot.",
        "`runZinbwave` is recommended by default.",
        sep = "\n"
    )
    weights <- tryCatch(
        expr = weights(object),
        error = function(e) {
            stop(msg, call. = FALSE)
        }
    )
    if (is.null(weights)) {
        stop(msg, call. = FALSE)
    }
    assert(is.matrix(weights))
    weights
}

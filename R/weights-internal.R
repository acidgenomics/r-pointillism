## Require zero weights to be calculated, showing a more informative error.
## Updated 2019-07-31.
.weights <- function(object) {
    msg <- paste(
        "Object does not contain zero weights in `weights` slot.",
        "`runZinbwave` is recommended by default.",
        sep = "\n"
    )
    weights <- tryCatch(
        expr = weights(object),
        error = function(e) {
            stop(msg)
        }
    )
    if (is.null(weights)) {
        stop(msg)
    }
    assert(is.matrix(weights))
    weights
}

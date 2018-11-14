# FIXME Improve appearance.
# Remove separator bar and import code from basejump.
# FIXME Switch to using `showSlotInfo()` (see bcbioRNASeq.)



#' @inherit methods::show
#' @name show
#' @author Michael Steinbuagh
NULL



#' @importFrom methods show
#' @aliases NULL
#' @export
methods::show



show.CellCycleMarkers <-  # nolint
    function(object) {
        validObject(object)

        # Include the organism information.
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["ensemblRelease"]]

        # Include the gene lengths per phase.
        split <- split(x = object, f = object[["phase"]])
        assert_that(is(split, "SplitDataFrameList"))
        lengths <- nrow(split)
        # FIXME Can we avoid `sapply()` here?
        genes <- sapply(
            X = split,
            FUN = function(x) {
                x <- x[["geneName"]]
                x <- sort(x)
                x <- c(head(x, n = 2L), "...", tail(x, n = 2L))
                paste(x, collapse = " ")
            },
            USE.NAMES = TRUE
        )

        return <- c(
            paste(class(object), metadata(object)[["version"]]),
            paste0(organism, " (Ensembl ", release, ")"),
            url,
            citation,
            separatorBar,
            paste0(names(genes), " (", lengths, "): ", genes)
        )

        cat(return, sep = "\n")
    }



show.CellTypeMarkers <-  # nolint
    function(object) {
        validObject(object)

        # Include the organism information.
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["ensemblRelease"]]

        # Include the gene lengths per phase.
        split <- split(x = object, f = object[["cellType"]])
        assert_that(is(split, "SplitDataFrameList"))
        lengths <- nrow(split)
        genes <- sapply(
            X = split,
            FUN = function(x) {
                x <- x[["geneName"]]
                x <- sort(x)
                x <- c(head(x, n = 2L), "...", tail(x, n = 2L))
                paste(x, collapse = " ")
            },
            USE.NAMES = TRUE
        )

        return <- c(
            paste(class(object), metadata(object)[["version"]]),
            paste0(organism, " (Ensembl ", release, ")"),
            url,
            citation,
            separatorBar,
            paste0(names(genes), " (", lengths, "): ", genes)
        )

        cat(return, sep = "\n")
    }



# FIXME Need to finish this and add S4 method.
show.Markers <-  # nolint
    function(object) {
        validObject(object)
        data <- as(object, "DataFrame")

        version <- metadata(object)[["version"]]
        date <- metadata(object)[["date"]]

        gr <- metadata(object)[["GRanges"]]
        organism <- metadata(gr)[["organism"]]
        build <- metadata(gr)[["build"]]
        release <- metadata(gr)[["release"]]

        if ("cluster" %in% colnames(data)) {
            f <- "FindAllMarkers"
            clusters <- levels(data[["cluster"]])
            split <- split(x = data, f = object[["cluster"]])
            assert_that(is(split, "SplitDataFrameList"))
            assert_are_identical(names(split), clusters)

        } else {
            f <- "FindMarkers"
        }

        stop("Not added yet.")
    }



#' #' @rdname show
#' #' @export
#' setMethod(
#'     f = "show",
#'     signature = signature("CellCycleMarkers"),
#'     definition = show.CellCycleMarkers
#' )




#' #' @rdname show
#' #' @export
#' setMethod(
#'     f = "show",
#'     signature = signature("CellTypeMarkers"),
#'     definition = show.CellTypeMarkers
#' )

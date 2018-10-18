globalVariables(".")

# TODO Make a header function in show.R
url <- "https://steinbaugh.com/pointillism"
citation <- "citation(\"pointillism\")"

# FIXME Use the function in call, it's more dynamic.
separatorBar <- basejump::separator()

# DR marker default color palettes
darkMarkerColors <-
    ggplot2::scale_color_viridis_c(option = "plasma")
lightMarkerColors <-
    ggplot2::scale_color_gradient(low = "gray90", high = "red")

continuousColor <- quote(getOption(
    "pointillism.continuous.color",
    ggplot2::scale_color_gradient(
        low = "gray75",
        high = "purple"
    )
))
continuousColor2 <- quote(getOption(
    "pointillism.continuous.color",
    ggplot2::scale_color_gradient(
        low = "orange",
        high = "purple"
    )
))
discreteColor <- quote(getOption("pointillism.discrete.color", NULL))

dark <- quote(getOption("pointillism.dark", FALSE))
# up, down, both (bcbioRNASeq).
dimsUse <- quote(c(1L, 2L))
direction <- c("positive", "negative", "both")
expression <- c("mean", "sum")
headerLevel <- 2L
label <- quote(getOption("pointillism.label", TRUE))
labelSize <- quote(getOption("pointillism.labelSize", 6L))
legend <- quote(getOption("pointillism.legend", TRUE))
pointAlpha <- quote(getOption("pointillism.pointAlpha", 0.8))
pointSize <- quote(getOption("pointillism.pointSize", 0.75))
pointsAsNumbers <- quote(getOption("pointillism.pointsAsNumbers", FALSE))
reducedDim <- "TSNE"

# Set default ggplot2 theme.
if (isTRUE(getOption("pointillism.dark"))) {
    theme <- basejump::theme_midnight
} else {
    theme <- basejump::theme_paperwhite
}
theme <- theme(base_size = 14L, legend_position = "right")
ggplot2::theme_set(new = theme)

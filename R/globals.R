globalVariables(".")

url <- "https://steinbaugh.com/pointillism"
citation <- "citation(\"pointillism\")"

separatorBar <- basejump::separator()

# DR marker default color palettes.
darkMarkerColors <-
    ggplot2::scale_colour_viridis_c(option = "plasma")
lightMarkerColors <-
    ggplot2::scale_colour_gradient(low = "gray90", high = "red")

# Recommend serial by default.
# This works reliably across platforms, but is slower.
BPPARAM <- quote(BiocParallel::SerialParam())  # nolint

continuousColor <- quote(getOption(
    "pointillism.continuous.color",
    ggplot2::scale_colour_gradient(
        low = "gray75",
        high = "purple"
    )
))
continuousColorPurpleOrange <- quote(getOption(
    "pointillism.continuous.color",
    ggplot2::scale_colour_gradient2(
        low = "orange",
        mid = "gray75",
        high = "purple",
        midpoint = 0L
    )
))

discreteColor <- quote(getOption("pointillism.discrete.color", default = NULL))

dark <- quote(getOption("pointillism.dark", default = FALSE))
# up, down, both (bcbioRNASeq).
dimsUse <- quote(c(1L, 2L))
direction <- c("up", "down", "both")
expression <- c("mean", "sum")
headerLevel <- 2L
label <- quote(getOption("pointillism.label", default = TRUE))
labelSize <- quote(getOption("pointillism.labelSize", default = 6L))
legend <- quote(getOption("pointillism.legend", default = TRUE))
pointAlpha <- quote(getOption("pointillism.pointAlpha", default = 0.85))
pointSize <- quote(getOption("pointillism.pointSize", default = 0.75))
pointsAsNumbers <-
    quote(getOption("pointillism.pointsAsNumbers", default = FALSE))
reducedDim <- "TSNE"

# Set default ggplot2 theme.
if (isTRUE(getOption("pointillism.dark"))) {
    theme <- minimalism::theme_midnight
} else {
    theme <- minimalism::theme_paperwhite
}
theme <- theme(base_size = 14L, legend_position = "right")
ggplot2::theme_set(new = theme)

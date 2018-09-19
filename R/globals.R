globalVariables(".")

separatorBar <- basejump::separator()
url <- "https://steinbaugh.com/pointillism"
citation <- "citation(\"pointillism\")"

# DR marker default color palettes
darkMarkerColors <-
    ggplot2::scale_color_viridis_c(option = "plasma")
lightMarkerColors <-
    ggplot2::scale_color_gradient(low = "gray90", high = "red")

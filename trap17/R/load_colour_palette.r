#' Load colour palette
#'
#' Load a custom, colour-blind-friendly palette of colours for results plotting
#' @export

load_colour_palette <- function() {

    res <- list()
    
    res[[1]] <- c("grey75",
                "#2f1855",
                "#98c339",
                "#6349bc",
                "#3dc662",
                "#5b006f",
                "#02da89",
                "#c759c8",
                "#32901a",
                "#e860c6",
                "#9aa400",
                "#0046a5",
                "#ffc746",
                "#43b3ff",
                "#cf5f18",
                "#02c0ad",
                "#8e0059",
                "#7fe9a5",
                "#78005f",
                "#bbe16f",
                "#60002f",
                "#74e8c5",
                "#740011",
                "#007f42",
                "#ff68a5",
                "#b4c87d",
                "#6e0027",
                "#ffbe58",
                "#e6b4ff",
                "#746600",
                "#ff7081",
                "#7f3d00")

    res[[2]] <- c("grey75",
                "#e79fff",
                "#ffad6a",
                "#b10b80",
                "#f2d250",
                "#d0b3ff",
                "#9669a9",
                "#015dcf",
                "#98000a",
                "#246a00",
                "#ff71ba",
                "#cb57c4",
                "#c1232e",
                "#d1dc70",
                "#59007b",
                "#ff6686",
                "#770055",
                "#fe8fff",
                "#01a0f2",
                "#5b50c4",
                "#391255",
                "#a99500",
                "#f34890",
                "#ff9276",
                "#ffafe6",
                "#d52a4c",
                "#cc4525",
                "#00b262",
                "#0160a6",
                "#ff9ec3",
                "#ffab4b",
                "#7c92d9")

    res[[3]] <- cbind(c("Clo", "Be", "Cap", "Cau", "En"),
                      c("#7c92d9", "#ffab4b", "#ff9ec3", "#0160a6", "#00b262"))

    return(res)
            
}

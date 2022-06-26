#' @title Shiny demo
#' @description Run a Shiny app which demonstrates the Kochanek-Bartels spline.
#'
#' @return No value returned.
#' @export
#' @importFrom shiny shinyAppDir
shinyKBS <- function(){
  appDir <- system.file("shiny", "threejs", package = "qsplines")
  shinyAppDir(appDir, options = list(launch.browser = TRUE))
}
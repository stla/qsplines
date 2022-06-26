library(shiny)

library(rgl)
options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = TRUE)

unitSphere <- subdivision3d(icosahedron3d(), depth = 1)
unitSphere$vb[4, ] <-
  apply(unitSphere$vb[1:3, ], 2, function(x) sqrt(sum(x * x)))
unitSphere$normals <- unitSphere$vb

drawSphere <- function(center=c(0, 0, 0), radius = 1, ...){
  sphere <- scale3d(unitSphere, radius, radius, radius)
  shade3d(translate3d(sphere, center[1], center[2], center[3]), ...)
}


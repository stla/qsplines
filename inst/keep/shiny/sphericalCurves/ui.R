library(shinycssloaders)

shinyUI(
  fluidPage(
    # titlePanel("Kochanek-Bartels spline"),
    sidebarLayout(
      sidebarPanel(
        numericInput("numt", "tension", value = -1, step = 0.1),
        numericInput("numc", "continuity", value = 5, step = 0.1),
        numericInput("numb", "bias", value = 0, step = 0.1),
        actionButton("test", "Test"),
      ),
      mainPanel(
        withSpinner(
          rglwidgetOutput("sphere")
        ),
      )
    )
  )
)
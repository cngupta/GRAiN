library(shiny)

source('global.R')


source('UI.R', local = TRUE)
source('server.R')

shinyApp(
  ui = UI,
  server = server
)
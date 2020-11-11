#' Open Shiny App
#'
#' Open the application graphical user interface to access `plink`
#' command and run easy statistical tests on your data.
#'
#' @export

run <- function() shiny::shinyApp(ui = ui, server = server)

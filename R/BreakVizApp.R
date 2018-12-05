library(shiny)
library(chromoMap)
library(shinythemes)
library(shinyalert)
# Define UI for application
ui <- fluidPage(
  #Set UI theme
  theme = shinytheme("spacelab"),
  #Set for using shinyalert
  useShinyalert(),
  #Application header
  headerPanel(
    h1(
      "Welcome to BreakViz",
      style = "font-family: 'Courier New', cursive;
      font-weight: 800; line-height: 1.1;
      color: #1905AF;"
    )
  ),
  #UI for sidebar
  sidebarLayout(
    #UI for choosing bedfile
    sidebarPanel(
      # UI for uploading input file
      fileInput(
        "bedFile",
        label = h3("Choose BED file"),
        multiple = FALSE,
        accept = c("text/comma-separated-values,text/plain", ".bed")
      ),
      # UI for changing color set
      selectInput(
        "colorSet",
        label = h3("Select color set"),
        choices = list(
          "green/blue" = 1,
          "red/yellow" = 2,
          "blue/yellow" = 3
        ),
        selected = 1
      ),
      #UI for sliders that change MinOverlapping value and MaxDistance value
      sliderInput(
        "MinOverlap",
        label = h3("Minimum basepairs that overlapped to two chromosomes (bp)"),
        min = 0,
        max = 1000,
        step = 10,
        value = 100
      ),
      sliderInput(
        "MaxDistance",
        label = h3("Maximum distance between two parts of one pair (bp)"),
        min = 0,
        max = 5000,
        value = 1000,
        step = 50
      )
    ),
    # Show a plot of the generated distribution
    mainPanel(chromoMapOutput("myChromoMap"))
  )
)

# Define server logic required to draw chromoMap
server <- function(input, output, session) {
  #download bedfile
  output$myChromoMap <- renderChromoMap({
    inputBedfile <- input$bedFile
    if (is.null(inputBedfile))
      return(NULL)
    #change bedfile format
    inputBedfile <-
      rtracklayer::import(inputBedfile$datapath, format = "bed")
    #catch error
    x <- tryCatch(
      visPossiblePair(
        inputBedfile,
        input$MinOverlap,
        input$MaxDistance,
        input$colorSet
      )
      ,
      error = function(e) {
        e
      }
    )
  })
  #update output by changing MinOverlap, MaxDistance and bedfile
  observeEvent({
    input$bedFile
    input$MaxDistance
    input$MinOverlap
  }, {
    inputBedfile <- input$bedFile
    if (is.null(inputBedfile))
      return(NULL)
    inputBedfile <-
      rtracklayer::import(inputBedfile$datapath, format = "bed")
    #catch error
    x <- tryCatch({
      visPossiblePair(inputBedfile,
                      input$MinOverlap,
                      input$MaxDistance,
                      input$colorSet)
    },
    error = function(e) {
      e
    })
    #collect error message
    mess <- x$message
    mess <- toString(mess)
    #send shinyalert motification
    if (mess != '') {
      shinyalert("Oops! Please change a new dataset", mess, type = "error")
    }
  })
}
# Run the application
shinyApp(ui = ui, server = server)

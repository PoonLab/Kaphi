library(shiny)

ui <- fluidPage(
  
  # Page Title
  titlePanel(
    h1(strong("Kaphi"), "- Kernel-embedded ABC-SMC for phylodynamic inference"),
    windowTitle = "Kaphi - Kernel-embedded ABC-SMC for phylodynamic inference"
  ),
  
  sidebarLayout(
    
    sidebarPanel( 
      # Row for Newick Text/File Input 
      fluidRow(
        h3(strong(em("Newick Input"))),
        textInput(inputId = "newickString", label = "Enter a Newick String"), 
        fileInput(inputId = "newickFile", label = "Choose a Newick File"),
        actionButton(inputId = "processString", label = "Process String"),
        actionButton(inputId = "processFile", label = "Process File")
      ),
      # Row for Configuration Creation
      fluidRow(
        h3(strong(em("Configuration Creation")))
      ),
      # Row for Model Specification
      fluidRow(
        h3(strong(em("Model Specification")))
      ),
      # Row for Running Simulation
      fluidRow(
        h3(strong(em("Running Simulation")))
      )
    ),
    
    mainPanel(
      tabsetPanel(
        # Tab for Tree Plot
        tabPanel("Tree Plot", plotOutput(outputId = "tree")), 
        # Tab for Prior Distributions
        tabPanel("Prior Distributions"), 
        # Tab for Feedback/Diagnosis
        tabPanel("Feedback/Diagnosis"),
        # Tab for Simulation Results
        tabPanel("Simulation Results")
      )
    )
    
  )
  
)  

server <- function(input, output) {
  
  newickInput <- reactiveValues(data = NULL)
  
  # Reading Tree from Newick String
  observeEvent(
    input$processString,
    {
      newickInput$data <- read.tree(text = input$newickString)
    }
  )
  
  # Reading Tree from Newick File
  observeEvent(
    input$processFile,
    {
      inFile <- input$newickFile
      newickInput$data <- read.tree(inFile$datapath)
    }
  )
  
  # Rendering Newick Input
  output$tree <- renderPlot({
    if (is.null(newickInput$data)) return()
    plot(newickInput$data)
  })
  
}

shinyApp(ui = ui, server = server)

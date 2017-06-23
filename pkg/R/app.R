library(shiny)
library(Kaphi)
library(rmarkdown)

ui <- fluidPage(
  
  # Page Title
  titlePanel(
    title = h1(strong("Kaphi"), "- Kernel-embedded ABC-SMC for phylodynamic inference"),
    windowTitle = "Kaphi - Kernel-embedded ABC-SMC for phylodynamic inference"
  ),
  
  sidebarLayout(
    
    sidebarPanel( 
      # Allowing Independent Scrolling in the Sidebar
      id = "sidebarPanel",
      style = "overflow-y:scroll; max-height:600px",
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
        h3(strong(em("SMC Settings Initialization"))),
        numericInput(inputId = "particleNumber", label = "Number of Particles", value = 1000),
        numericInput(inputId = "sampleNumber", label = "Number of Samples", value = 5),
        numericInput(inputId = "ESSTolerance", label = "Effective Sample Size (ESS) Tolerance", value = 1.5),
        numericInput(inputId = "finalEpsilon", label = "Final Epsilon", value = 0.01),
        numericInput(inputId = "finalAcceptanceRate", label = "Final Acceptance Rate", value = 0.015),
        numericInput(inputId = "quality", label = "Quality", value = 0.95),
        numericInput(inputId = "stepTolerance", label = "Step Tolerance", value = 1e-5),
        actionButton(inputId = "initializeSMCSettings", label = "Initialize SMC Settings")
      ),
      # Row for Model Selection and Initialization
      fluidRow(
        h3(strong(em("Model Selection and Initialization")))
      ),
      # Row for Running Simulation
      fluidRow(
        h3(strong(em("Running Simulation")))
      )
    ),
    
    mainPanel(
      tabsetPanel(
        # Tab for Tree Plot
        tabPanel(
          title = "Tree Plot",
          textInput(inputId = "treeTitle", label = "Enter Tree Title"),
          fluidRow(
            column(
              6,
              sliderInput("width", "Plot Width (px)", min = 0, max = 10000, value = 500)
            ),
            column(
              6,
              sliderInput("height", "Plot Height (px)", min = 0, max = 10000, value = 500)
            )
          ),
          uiOutput("tree.ui")
        ), 
        # Tab for Prior Distributions
        tabPanel(title = "Prior Distributions"), 
        # Tab for Feedback/Diagnosis
        tabPanel(title = "Feedback/Diagnosis"),
        # Tab for Simulation Results
        tabPanel(title = "Simulation Results")
      )
    )
    
  )
  
)  

server <- function(input, output) {
  
  newickInput <- reactiveValues(data = NULL)
  
  config <- list(
    params=NA,
    priors=list(),
    prior.densities=list(),
    proposals=list(),
    proposal.densities=list(),
    model=NA,
    
    # SMC settings
    nparticle=1000,
    nsample=5,
    ess.tolerance=1.5,
    final.epsilon=0.01,
    final.accept.rate=0.015,
    quality=0.95,
    step.tolerance=1e-5,
    
    # kernel settings
    decay.factor=0.2,
    rbf.variance=2.0,
    sst.control=1.0,
    norm.mode='MEAN'
  )
  
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
  
  # Plotting Newick Input
  output$tree <- renderPlot({
    if (is.null(newickInput$data)) return()
    plot(newickInput$data, main = input$treeTitle)
  })
  
  # Rendering Newick Input
  output$tree.ui <- renderUI({
    plotOutput("tree", width = input$width, height = input$height)
  })
  
  # Initializing SMC Settings
  observeEvent(
    input$initializeSMCSettings,
    {
      config$nparticle <- input$particleNumber
      config$nsample <- input$sampleNumber
      config$ess.tolerance <- input$ESSTolerance
      config$final.epsilon <- input$finalEpsilon
      config$final.accept.rate <- input$finalAcceptanceRate
      config$quality <- input$quality
      config$step.tolerance <- input$stepTolerance
    }
  )
  
}

shinyApp(ui = ui, server = server)

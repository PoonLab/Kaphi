library(shiny)

ui <- fluidPage(
  
  # Page Title
  titlePanel(
    title = h1(strong("Kaphi"), "- Kernel-embedded ABC-SMC for phylodynamic inference"),
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
        h3(strong(em("Model Specification"))),
        # General Model Drop-down Menu
        selectInput(
          inputId = "generalModel", 
          label = "Select General Simulation Model", 
          choices = c(
            Coalescent = "coalescent",
            Compartmental = "compartmental", 
            Network = "network",
            Speciation = "speciation"
          )
        ),
        # Specific Model Drop-down Menus
        # Coalescent Model Drop-down Menu
        conditionalPanel(
          condition = "input.generalModel == 'coalescent'",
          selectInput(
            inputId = "specificCoalescent", 
            label = "Select Specific Coalescent Model", 
            choices = c(
              "Constant Coalescent"
            )
          )
        ),
        # Coalescent Model Drop-down Menu
        conditionalPanel(
          condition = "input.generalModel == 'compartmental'",
          selectInput(
            inputId = "specificCompartmental", 
            label = "Select Specific Compartmental Model", 
            choices = c(
              "Susceptible-Infected-Removed-Dynamic (SIRD)",
              "Susceptible-Infected-Removed-Non-Dynamic (SIRND)",
              "Susceptible-Exposed-Infected-Removed (SEIR)",
              "Susceptible-Infected-Susceptible (SIS)"
            )
          )
        ),
        # Coalescent Model Drop-down Menu
        conditionalPanel(
          condition = "input.generalModel == 'network'",
          selectInput(
            inputId = "specificNetwork", 
            label = "Select Specific Network Model", 
            choices = c(
              "Coming Soon"
            )
          )
        ),
        # Coalescent Model Drop-down Menu
        conditionalPanel(
          condition = "input.generalModel == 'speciation'",
          selectInput(
            inputId = "specificSpeciation", 
            label = "Select Specific Speciation Model", 
            choices = c(
              "Yule", 
              "Birth-Death",
              "Binary State Speciation Extinction (BiSSE)",
              "MuSSE",
              "QuaSSE",
              "GeoSSE",
              "BiSS-ness",
              "ClaSSE"
            )
          )
        )
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
          selectInput(inputId = "downloadFormat", label = "Select Download Format", choices = c(PNG = "png", PDF = "pdf")),
          downloadButton(outputId = "downloadTree", label = "Download Tree"),
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
  
  # Downloading Tree Plot
  output$downloadTree <- downloadHandler(
    fileName <-  function() {
      paste(input$treeTitle, input$downloadFormat, sep=".")
    },
    content <- function(file) {
      if(input$downloadFormat == "png")png(file) # open the png device
      else pdf(file) # open the pdf device
      plot(newickInput$data, main = input$treeTitle) # draw the plot
      dev.off()  # turn the device off
    }
  )
  
}

shinyApp(ui = ui, server = server)

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
        h3(strong(em("SMC Settings Initialization"))),
        numericInput(inputId = "particleNumber", label = "Number of Particles", value = 1000),
        numericInput(inputId = "sampleNumber", label = "Number of Samples", value = 5),
        numericInput(inputId = "ESSTolerance", label = "Effective Sample Size (ESS) Tolerance", value = 1.5),
        numericInput(inputId = "finalEpsilon", label = "Final Epsilon", value = 0.01),
        numericInput(inputId = "finalAcceptanceRate", label = "Final Acceptance Rate", value = 0.015),
        numericInput(inputId = "quality", label = "Quality", value = 0.95),
        numericInput(inputId = "stepTolerance", label = "Step Tolerance", value = 1e-5)
      ),
      # Row for Model Selection and Initialization
      fluidRow(
        h3(strong(em("Model Selection and Initialization"))),
        # General Model Drop-down Menu
        selectInput(
          inputId = "generalModel", 
          label = "Select General Simulation Model", 
          choices = c(
            "Coalescent",
            "Compartmental", 
            "Network",
            "Speciation"
          )
        ),
        # Specific Model Drop-down Menus
        # Coalescent Model Drop-down Menu
        conditionalPanel(
          condition = "input.generalModel == 'Coalescent'",
          selectInput(
            inputId = "specificCoalescent", 
            label = "Select Specific Coalescent Model", 
            choices = c(
              "Constant Coalescent"
            )
          ),
          conditionalPanel(
            condition = "input.specificCoalescent == 'Constant Coalescent'",
            selectInput(
              inputId = "NeTauPriorDistribution", 
              label = "Ne Tau Prior Distribution",  
              choices = c(
                "Exponential",
                "Gamma",
                "Normal",
                "Log Normal"
              )
            ),
            # conditionalPanel(
            #   condition = "input.NeTauPriorDistribution == 'Exponential'",
            #   numericInput(inputId = "NeTauPriorExponentialRate", label = "Ne Tau Prior Exponential Rate", value = )
            # ),
            # conditionalPanel(
            #   condition = "input.NeTauPriorDistribution == 'Gamma'",
            #   numericInput(inputId = "NeTauPriorGammaShape", label = "Ne Tau Prior Gamma Shape", value = ),
            #   numericInput(inputId = "NeTauPriorGammaRate", label = "Ne Tau Prior Gamma Rate", value = )
            # ),
            # conditionalPanel(
            #   condition = "input.NeTauPriorDistribution == 'Normal'",
            #   numericInput(inputId = "NeTauPriorNormalMean", label = "Ne Tau Prior Normal Mean", value = ),
            #   numericInput(inputId = "NeTauPriorNormalStandardDeviation", label = "Ne Tau Prior Normal Standard Deviation", value = )
            # ),
            # conditionalPanel(
            #   condition = "input.NeTauPriorDistribution == 'Log Normal'",
            #   numericInput(inputId = "NeTauPriorLogNormalMean", label = "Ne Tau Prior Log Normal Mean", value = ),
            #   numericInput(inputId = "NeTauPriorLogNormalStandardDeviation", label = "Ne Tau Prior Log Normal Standard Deviation", value = )
            # ),
            selectInput(
              inputId = "NeTauProposalDistribution", 
              label = "Ne Tau Proposal Distribution", 
              choices = c(
                "Exponential",
                "Gamma",
                "Normal",
                "Log Normal"
              )
              # ,
              # conditionalPanel(
              #   condition = "input.NeTauProposalDistribution == 'Exponential'",
              #   numericInput(inputId = "NeTauProposalExponentialRate", label = "Ne Tau Proposal Exponential Rate", value = )
              # ),
              # conditionalPanel(
              #   condition = "input.NeTauProposalDistribution == 'Gamma'",
              #   numericInput(inputId = "NeTauProposalGammaShape", label = "Ne Tau Proposal Gamma Shape", value = ),
              #   numericInput(inputId = "NeTauProposalGammaRate", label = "Ne Tau Proposal Gamma Rate", value = )
              # ),
              # conditionalPanel(
              #   condition = "input.NeTauProposalDistribution == 'Normal'",
              #   numericInput(inputId = "NeTauProposalNormalMean", label = "Ne Tau Proposal Normal Mean", value = ),
              #   numericInput(inputId = "NeTauProposalNormalStandardDeviation", label = "Ne Tau Proposal Normal Standard Deviation", value = )
              # ),
              # conditionalPanel(
              #   condition = "input.NeTauProposalDistribution == 'Log Normal'",
              #   numericInput(inputId = "NeTauProposalLogNormalMean", label = "Ne Tau Proposal Log Normal Mean", value = ),
              #   numericInput(inputId = "NeTauProposalLogNormalStandardDeviation", label = "Ne Tau Proposal Log Normal Standard Deviation", value = )
              # )
            )
          )
        ),
        # Compartmental Model Drop-down Menu
        conditionalPanel(
          condition = "input.generalModel == 'Compartmental'",
          selectInput(
            inputId = "specificCompartmental", 
            label = "Select Specific Compartmental Model", 
            choices = c(
              "Susceptible-Infected-Removed-Dynamic (SIRD)",
              "Susceptible-Infected-Removed-Non-Dynamic (SIRND)",
              "Susceptible-Exposed-Infected-Removed (SEIR)",
              "Susceptible-Infected-Susceptible (SIS)"
            )
          ),
          conditionalPanel(
            condition = "input.specificCompartmental == 'Susceptible-Infected-Removed-Dynamic (SIRD)'"
          ),
          conditionalPanel(
            condition = "input.specificCompartmental == 'Susceptible-Infected-Removed-Non-Dynamic (SIRND)'"
          ),
          conditionalPanel(
            condition = "input.specificCompartmental == 'Susceptible-Exposed-Infected-Removed (SEIR)'"
          ),
          conditionalPanel(
            condition = "input.specificCompartmental == 'Susceptible-Infected-Susceptible (SIS)'"
          )
        ),
        # Network Model Drop-down Menu
        conditionalPanel(
          condition = "input.generalModel == 'Network'",
          selectInput(
            inputId = "specificNetwork",
            label = "Select Specific Network Model",
            choices = c(
              "Coming Soon"
            )
          ),
          conditionalPanel(
            condition = "input.specificNetwork == 'Coming Soon'"
          )
        ),
        # Speciation Model Drop-down Menu
        conditionalPanel(
          condition = "input.generalModel == 'Speciation'",
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
          ),
          conditionalPanel(
            condition = "input.specificSpeciation == 'Yule'",
            selectInput(
              inputId = "lambdaPriorDistribution", 
              label = "Lambda Prior Distribution",  
              choices = c(
                "Exponential",
                "Gamma",
                "Normal",
                "Log Normal"
              )
            ),
            # conditionalPanel(
            #   condition = "input.lambdaPriorDistribution == 'Exponential'",
            #   numericInput(inputId = , label = , value = )
            # ),
            # conditionalPanel(
            #   condition = "input.lambdaPriorDistribution == 'Gamma'",
            #   numericInput(inputId = , label = , value = ),
            #   numericInput(inputId = , label = , value = )
            # ),
            # conditionalPanel(
            #   condition = "input.lambdaPriorDistribution == 'Normal'",
            #   numericInput(inputId = , label = , value = ),
            #   numericInput(inputId = , label = , value = )
            # ),
            # conditionalPanel(
            #   condition = "input.lambdaPriorDistribution == 'Log Normal'",
            #   numericInput(inputId = , label = , value = ),
            #   numericInput(inputId = , label = , value = )
            # ),
            selectInput(
              inputId = "lambdaProposalDistribution", 
              label = "Lambda Proposal Distribution", 
              choices = c(
                "Exponential",
                "Gamma",
                "Normal",
                "Log Normal"
              )
            )
            # ,
            # conditionalPanel(
            #   condition = "input.lambdaProposalDistribution == 'Exponential'",
            #   numericInput(inputId = , label = , value = )
            # ),
            # conditionalPanel(
            #   condition = "input.lambdaProposalDistribution == 'Gamma'",
            #   numericInput(inputId = , label = , value = ),
            #   numericInput(inputId = , label = , value = )
            # ),
            # conditionalPanel(
            #   condition = "input.lambdaProposalDistribution == 'Normal'",
            #   numericInput(inputId = , label = , value = ),
            #   numericInput(inputId = , label = , value = )
            # ),
            # conditionalPanel(
            #   condition = "input.lambdaProposalDistribution == 'Log Normal'",
            #   numericInput(inputId = , label = , value = ),
            #   numericInput(inputId = , label = , value = )
            # )
          ),
          conditionalPanel(
            condition = "input.specificSpeciation == 'Birth-Death'"
          ),
          conditionalPanel(
            condition = "input.specificSpeciation == 'Binary State Speciation Extinction (BiSSE)'"
          ),
          conditionalPanel(
            condition = "input.specificSpeciation == 'MuSSE'"
          ),
          conditionalPanel(
            condition = "input.specificSpeciation == 'QuaSSE'"
          ),
          conditionalPanel(
            condition = "input.specificSpeciation == 'GeoSSE'"
          ),
          conditionalPanel(
            condition = "input.specificSpeciation == 'BiSS-ness'"
          ),
          conditionalPanel(
            condition = "input.specificSpeciation == 'ClaSSE'"
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

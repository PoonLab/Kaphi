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
        numericInput(inputId = "stepTolerance", label = "Step Tolerance", value = 1e-5),
        actionButton(inputId = "initializeSMCSettings", label = "Initialize SMC Settings")
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
            conditionalPanel(
              condition = "input.NeTauPriorDistribution == 'Exponential'",
              numericInput(inputId = "NeTauCoalescentPriorExponentialRate", label = "Ne Tau Prior Exponential Rate", value = 0),
              actionButton(inputId = "initializeNeTauCoalescentPriorExponential", label = "Initialize Ne Tau Prior")
            ),
            conditionalPanel(
              condition = "input.NeTauPriorDistribution == 'Gamma'",
              numericInput(inputId = "NeTauCoalescentPriorGammaShape", label = "Ne Tau Prior Gamma Shape", value = 0),
              numericInput(inputId = "NeTauCoalescentPriorGammaRate", label = "Ne Tau Prior Gamma Rate", value = 0),
              actionButton(inputId = "initializeNeTauCoalescentPriorGamma", label = "Initialize Ne Tau Prior")
            ),
            conditionalPanel(
              condition = "input.NeTauPriorDistribution == 'Normal'",
              numericInput(inputId = "NeTauCoalescentPriorNormalMean", label = "Ne Tau Prior Normal Mean", value = 0),
              numericInput(inputId = "NeTauCoalescentPriorNormalStandardDeviation", label = "Ne Tau Prior Normal Standard Deviation", value = 0),
              actionButton(inputId = "initializeNeTauCoalescentPriorNormal", label = "Initialize Ne Tau Prior")
            ),
            conditionalPanel(
              condition = "input.NeTauPriorDistribution == 'Log Normal'",
              numericInput(inputId = "NeTauCoalescentPriorLogNormalMean", label = "Ne Tau Prior Log Normal Mean", value = 0),
              numericInput(inputId = "NeTauCoalescentPriorLogNormalStandardDeviation", label = "Ne Tau Prior Log Normal Standard Deviation", value = 0),
              actionButton(inputId = "initializeNeTauCoalescentPriorLogNormal", label = "Initialize Ne Tau Prior")
            ),
            selectInput(
              inputId = "NeTauProposalDistribution", 
              label = "Ne Tau Proposal Distribution", 
              choices = c(
                "Exponential",
                "Gamma",
                "Normal",
                "Log Normal"
              )
            ),
            conditionalPanel(
              condition = "input.NeTauProposalDistribution == 'Exponential'",
              numericInput(inputId = "NeTauCoalescentProposalExponentialRate", label = "Ne Tau Proposal Exponential Rate", value = 0),
              actionButton(inputId = "initializeNeTauCoalescentProposalExponential", label = "Initialize Ne Tau Proposal")
            ),
            conditionalPanel(
              condition = "input.NeTauProposalDistribution == 'Gamma'",
              numericInput(inputId = "NeTauCoalescentProposalGammaShape", label = "Ne Tau Proposal Gamma Shape", value = 0),
              numericInput(inputId = "NeTauCoalescentProposalGammaRate", label = "Ne Tau Proposal Gamma Rate", value = 0),
              actionButton(inputId = "initializeNeTauCoalescentProposalGamma", label = "Initialize Ne Tau Proposal")
            ),
            conditionalPanel(
              condition = "input.NeTauProposalDistribution == 'Normal'",
              numericInput(inputId = "NeTauCoalescentProposalNormalMean", label = "Ne Tau Proposal Normal Mean", value = 0),
              numericInput(inputId = "NeTauCoalescentProposalNormalStandardDeviation", label = "Ne Tau Proposal Normal Standard Deviation", value = 0),
              actionButton(inputId = "initializeNeTauCoalescentProposalNormal", label = "Initialize Ne Tau Proposal")
            ),
            conditionalPanel(
              condition = "input.NeTauProposalDistribution == 'Log Normal'",
              numericInput(inputId = "NeTauCoalescentProposalLogNormalMean", label = "Ne Tau Proposal Log Normal Mean", value = 0),
              numericInput(inputId = "NeTauCoalescentProposalLogNormalStandardDeviation", label = "Ne Tau Proposal Log Normal Standard Deviation", value = 0),
              actionButton(inputId = "initializeNeTauCoalescentProposalLogNormal", label = "Initialize Ne Tau Proposal")
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
            conditionalPanel(
              condition = "input.lambdaPriorDistribution == 'Exponential'",
              numericInput(inputId = "lambdaYulePriorExponentialRate", label = "Lambda Prior Exponential Rate", value = 0),
              actionButton(inputId = "initializeLambdaYulePriorExponential", label = "Initialize Lambda Prior")
            ),
            conditionalPanel(
              condition = "input.lambdaPriorDistribution == 'Gamma'",
              numericInput(inputId = "lambdaYulePriorGammaShape", label = "Lambda Prior Gamma Shape", value = 0),
              numericInput(inputId = "lambdaYulePriorGammaRate", label = "Lambda Prior Gamma Rate", value = 0),
              actionButton(inputId = "initializeLambdaYulePriorGamma", label = "Initialize Lambda Prior")
            ),
            conditionalPanel(
              condition = "input.lambdaPriorDistribution == 'Normal'",
              numericInput(inputId = "lambdaYulePriorNormalMean", label = "Lambda Prior Normal Mean", value = 0),
              numericInput(inputId = "lambdaYulePriorNormalStandardDeviation", label = "Lambda Prior Normal Standard Deviation", value = 0),
              actionButton(inputId = "initializeLambdaYulePriorNormal", label = "Initialize Lambda Prior")
            ),
            conditionalPanel(
              condition = "input.lambdaPriorDistribution == 'Log Normal'",
              numericInput(inputId = "lambdaYulePriorLogNormalMean", label = "Lambda Prior Log Normal Mean", value = 0),
              numericInput(inputId = "lambdaYulePriorLogNormalStandardDeviation", label = "Lambda Prior Log Normal Standard Deviation", value = 0),
              actionButton(inputId = "initializeLambdaYulePriorNormal", label = "Initialize Lambda Prior")
            ),
            selectInput(
              inputId = "lambdaProposalDistribution", 
              label = "Lambda Proposal Distribution", 
              choices = c(
                "Exponential",
                "Gamma",
                "Normal",
                "Log Normal"
              )
            ),
            conditionalPanel(
              condition = "input.lambdaProposalDistribution == 'Exponential'",
              numericInput(inputId = "lambdaYuleProposalExponentialRate", label = "Lambda Proposal Exponential Rate", value = 0),
              actionButton(inputId = "initializeLambdaYuleProposalExponential", label = "Initialize Lambda Proposal")
            ),
            conditionalPanel(
              condition = "input.lambdaProposalDistribution == 'Gamma'",
              numericInput(inputId = "lambdaYuleProposalGammaShape", label = "Lambda Proposal Gamma Shape", value = 0),
              numericInput(inputId = "lambdaYuleProposalGammaRate", label = "Lambda Proposal Gamma Rate", value = 0),
              actionButton(inputId = "initializeLambdaYuleProposalGamma", label = "Initialize Lambda Proposal")
            ),
            conditionalPanel(
              condition = "input.lambdaProposalDistribution == 'Normal'",
              numericInput(inputId = "lambdaYuleProposalNormalMean", label = "Lambda Proposal Normal Mean", value = 0),
              numericInput(inputId = "lambdaYuleProposalNormalStandardDeviation", label = "Lambda Proposal Normal Standard Deviation", value = 0),
              actionButton(inputId = "initializeLambdaYuleProposalNormal", label = "Initialize Lambda Proposal")
            ),
            conditionalPanel(
              condition = "input.lambdaProposalDistribution == 'Log Normal'",
              numericInput(inputId = "lambdaYuleProposalLogNormalMean", label = "Lambda Proposal Log Normal Mean", value = 0),
              numericInput(inputId = "lambdaYuleProposalLogNormalStandardDeviation", label = "Lambda Proposal Log Normal Standard Deviation", value = 0),
              actionButton(inputId = "initializeLambdaYuleProposalLogNormal", label = "Initialize Lambda Proposal")
            )
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

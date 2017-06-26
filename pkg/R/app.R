library(shiny)

distributions = list(
  "exp" = list(
    "rate" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    )
  ),
  "gamma" = list(
    "rate" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    ),
    "shape" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    )
  ),
  "lnorm" = list(
    "mean" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    ),
    "sd" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    )
  ),
  "norm" = list(
    "mean" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    ),
    "sd" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    )
  )
)

ConstantCoalescent = list(
  Ne.tau = distributions
)

SIRD = list(
  beta = distributions,
  gamma = distributions,
  mu = distributions
)

SIRND =list(
  beta = distributions,
  gamma = distributions
)

SIS = list(
  beta = distributions,
  gamma = distributions,
  mu = distributions
)

SEIR = list(
  beta = distributions,
  gamma = distributions,
  mu = distributions, 
  alpha = distributions
)

Yule = list(
  lambda = distributions
)

BirthDeath = list(
  lambda = distributions,
  mu = distributions
)  

BiSSE = list(
  lambda0 = distributions,
  lambda1 = distributions,
  mu0 = distributions,
  mu1 = distributions,
  q01 = distributions,
  q10 = distributions
)

MuSSE = list(
  lambda1 = distributions,
  lambda2 = distributions,
  lambda3 = distributions,
  mu1 = distributions,
  mu2 = distributions,
  mu3 = distributions,
  q12 = distributions,
  q13 = distributions,
  q21 = distributions,
  q23 = distributions,
  q31 = distributions,
  q32 = distributions
)

QuaSSE = list(
  lambda = distributions,
  mu = distributions,
  char = distributions
) 

GeoSSE = list(
  sA = distributions,
  sB = distributions,
  sAB = distributions,
  xA = distributions,
  xB = distributions,
  dA = distributions,
  dB = distributions
)

BiSSness = list(
  lambda0 = distributions,
  lambda1 = distributions,
  mu0 = distributions,
  mu1 = distributions,
  q01 = distributions,
  q10 = distributions,
  p0c = distributions,
  p0a = distributions,
  p1c = distributions,
  p1a = distributions
)

ClaSSE = list(
  lambda111 = distributions,
  lambda112 = distributions,
  lambda122 = distributions,
  lambda211 = distributions,
  lambda212 = distributions,
  lambda222 = distributions,
  mu1 = distributions,
  mu2 = distributions,
  q12 = distributions,
  q21 = distributions
)

models = list(
  "Coalescent" = list(
    "Constant Coalescent" = list(
      "Priors" = ConstantCoalescent,
      "Proposals" = ConstantCoalescent
    )
  ),
  "Compartmental" = list(
    "Susceptible-Infected-Removed-Dynamic (SIRD)" = list(
      "Priors" = SIRD,
      "Proposals" = SIRD
    ),
    "Susceptible-Infected-Removed-Non-Dynamic (SIRND)" = list(
      "Priors" = SIRND,
      "Proposals" = SIRND
    ),
    "Susceptible-Exposed-Infected-Removed (SEIR)" = list(
      "Priors" = SEIR,
      "Proposals" = SEIR
    ),
    "Susceptible-Infected-Susceptible (SIS)" = list(
      "Priors" = SIS,
      "Proposals" = SIS
    )
  ),
  "Networks" = list(),
  "Speciation" = list(
    "Yule" = list(
      "Priors" = Yule,
      "Proposals" = Yule
    ), 
    "Birth-Death" = list(
      "Priors" = BirthDeath,
      "Proposals" = BirthDeath
    ),
    "Binary State Speciation Extinction (BiSSE)" = list(
      "Priors" = BiSSE,
      "Proposals" = BiSSE
    ),
    "MuSSE" = list(
      "Priors" = MuSSE,
      "Proposals" = MuSSE
    ),
    "QuaSSE" = list(
      "Priors" = QuaSSE,
      "Proposals" = QuaSSE
    ),
    "GeoSSE" = list(
      "Priors" = GeoSSE,
      "Proposals" = GeoSSE
    ),
    "BiSS-ness" = list(
      "Priors" = BiSSness,
      "Proposals" = BiSSness
    ),
    "ClaSSE" = list(
      "Priors" = ClaSSE,
      "Proposals" = ClaSSE
    )
  )
)

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

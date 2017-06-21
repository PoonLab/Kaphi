distributions = list(
  "Exponential" = list(
    "Rate" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    )
  ),
  "Gamma" = list(
    "Rate" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    ),
    "Shape" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    )
  ),
  "Log Normal" = list(
    "Log Mean" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    ),
    "Log Standard Deviation" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    )
  ),
  "Normal" = list(
    "Mean" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    ),
    "Standard Deviation" = list(
      "Lower" = 0,
      "Upper" = Inf,
      "Default" = 1
    )
  )
)

models = list(
  "Coalescent" = list(
    "Constant Coalescent" = list(
      "Priors" = list(),
      "Proposals" = list()
    )
  ),
  "Compartmental" = list(
    "Susceptible-Infected-Removed-Dynamic (SIRD)" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "Susceptible-Infected-Removed-Non-Dynamic (SIRND)" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "Susceptible-Exposed-Infected-Removed (SEIR)" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "Susceptible-Infected-Susceptible (SIS)" = list(
      "Priors" = list(),
      "Proposals" = list()
    )
  ),
  "Networks" = list(),
  "Speciation" = list(
    "Yule" = list(
      "Priors" = list(),
      "Proposals" = list()
    ), 
    "Birth-Death" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "Binary State Speciation Extinction (BiSSE)" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "MuSSE" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "QuaSSE" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "GeoSSE" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "BiSS-ness" = list(
      "Priors" = list(),
      "Proposals" = list()
    ),
    "ClaSSE" = list(
      "Priors" = list(),
      "Proposals" = list()
    )
  )
)


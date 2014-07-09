module Brownian

using Distributions
using PDMats

import Base: rand!, rand
import Distributions: VariateForm, Univariate, Multivariate, ValueSupport, Discrete, Continuous, Sampleable

export
  StochasticProcess,
  UnivariateStochasticProcess,
  MultivariateStochasticProcess,
  DiscreteStochasticProcess,
  ContinuousStochasticProcess,
  #DiscreteUnivariateStochasticProcess,
  ContinuousUnivariateStochasticProcess,
  #DiscreteMultivariateStochasticProcess,
  ContinuousMultivariateStochasticProcess,
  BrownianMotion,
  #MvBrownianMotion,
  FGN,
  FBM,
  cov,
  rand!,
  rand

include("StochasticProcess.jl")
include(joinpath("univariate", "BrownianMotion.jl"))
include(joinpath("univariate", "FBM.jl"))
#include(joinpath("multivariate", "MvBrownianMotion.jl"))

end # module

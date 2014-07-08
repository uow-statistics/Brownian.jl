module Brownian

using Distributions
using PDMats

import Base: rand!, rand
import Distributions: VariateForm, Univariate, Multivariate, ValueSupport, Discrete, Continuous

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
  MvBrownianMotion,
  rand!,
  rand

include(joinpath("processes", "StochasticProcess.jl"))
include(joinpath("processes", "univariate", "BrownianMotion.jl"))
include(joinpath("processes", "univariate", "FBM.jl"))
include(joinpath("processes", "multivariate", "MvBrownianMotion.jl"))

end # module

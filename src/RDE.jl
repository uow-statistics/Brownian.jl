module RDE

using Distributions

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
  #ContinuousMultivariateStochasticProcess,
  BrownianMotion,
  rand!,
  rand

include(joinpath("processes", "StochasticProcess.jl"))
include(joinpath("processes", "univariate", "BrownianMotion.jl"))

end # module

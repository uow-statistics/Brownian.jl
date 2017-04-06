__precompile__()

module Brownian

using Distributions
using PDMats
using StatsBase

import Base: convert, rand!, rand
import Distributions: VariateForm, Univariate, Multivariate, ValueSupport, Discrete, Continuous, Sampleable
import StatsBase: IntegerVector, autocov!, autocov

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
  FBM,
  FGN,
  convert,
  autocov!,
  autocov,
  chol_update,
  rand!,
  rand

include("StochasticProcess.jl")
include(joinpath("univariate", "BrownianMotion.jl"))
include(joinpath("univariate", "FBM.jl"))
#include(joinpath("multivariate", "MvBrownianMotion.jl"))

end # module

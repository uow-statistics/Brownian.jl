module Brownian

using Distributions
using PDMats
using StatsBase

import Base: rand!, rand
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
  autocov!,
  autocov,
  rand!,
  rand,
  rand_fft

include("StochasticProcess.jl")
include(joinpath("univariate", "BrownianMotion.jl"))
include(joinpath("univariate", "FBM.jl"))
#include(joinpath("multivariate", "MvBrownianMotion.jl"))

end # module

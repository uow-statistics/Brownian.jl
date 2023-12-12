__precompile__()

module Brownian

using Distributions
using PDMats
using StatsBase
using Random
using FFTW
using LinearAlgebra
using SpecialFunctions

import LinearAlgebra: cholesky
import SpecialFunctions: gamma
import FFTW: fft
import Base: convert, rand
import Random: rand!
import Distributions:
    VariateForm, Univariate, Multivariate, ValueSupport, Discrete, Continuous, Sampleable
import StatsBase: autocov!, autocov

export StochasticProcess,
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

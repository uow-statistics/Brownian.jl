abstract type StochasticProcess{F<:VariateForm,S<:ValueSupport} <: Sampleable{F,S} end

const UnivariateStochasticProcess{S<:ValueSupport} = StochasticProcess{Univariate,S}
const MultivariateStochasticProcess{S<:ValueSupport} = StochasticProcess{Multivariate,S}

const DiscreteStochasticProcess{F<:VariateForm} = StochasticProcess{F,Discrete}
const ContinuousStochasticProcess{F<:VariateForm} = StochasticProcess{F,Continuous}

#const DiscreteUnivariateStochasticProcess = StochasticProcess{Univariate, Discrete}
const ContinuousUnivariateStochasticProcess = StochasticProcess{Univariate,Continuous}
#const DiscreteMultivariateStochasticProcess = StochasticProcess{Multivariate, Discrete}
const ContinuousMultivariateStochasticProcess = StochasticProcess{Multivariate,Continuous}

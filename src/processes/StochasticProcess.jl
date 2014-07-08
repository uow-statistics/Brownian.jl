abstract StochasticProcess{F<:VariateForm, S<:ValueSupport} <: Sampleable{F,S}

typealias UnivariateStochasticProcess{S<:ValueSupport} StochasticProcess{Univariate, S}
typealias MultivariateStochasticProcess{S<:ValueSupport} StochasticProcess{Multivariate, S}

typealias DiscreteStochasticProcess{F<:VariateForm} StochasticProcess{F, Discrete}
typealias ContinuousStochasticProcess{F<:VariateForm} StochasticProcess{F, Continuous}

#typealias DiscreteUnivariateStochasticProcess StochasticProcess{Univariate, Discrete}
typealias ContinuousUnivariateStochasticProcess StochasticProcess{Univariate, Continuous}
#typealias DiscreteMultivariateStochasticProcess StochasticProcess{Multivariate, Discrete}
typealias ContinuousMultivariateStochasticProcess StochasticProcess{Multivariate, Continuous}

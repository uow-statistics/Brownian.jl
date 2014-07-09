immutable MvBrownianMotion{Cov<:AbstractPDMat} <: ContinuousMultivariateStochasticProcess
  timepoints::Vector{Float64}
  npoints::Int64
  Σ::Cov

  function MvBrownianMotion{Cov<:AbstractPDMat}(t::Vector{Float64}, n::Int64, Σ::Cov)
    t[1] > 0.0 || error("First provided time point must be positive for Brownian motion.")
    issorted(t, lt=<=) || error("The time points must be strictly sorted.")
    int64(length(t)) == n || error("Number of time points must be equal to the vector holding the time points.")
    MvBrownianMotion{Cov<:AbstractPDMat}(t, n, Σ)
  end
end

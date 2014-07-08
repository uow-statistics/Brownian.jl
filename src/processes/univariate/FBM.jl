immutable FBM <: ContinuousUnivariateStochasticProcess
  timepoints::Vector{Float64}
  npoints::Int64
  hurst::Float64

  function FBM(t::Vector{Float64}, n::Int64)
    t[1] > 0.0 || error("Second time point must be positive for Brownian motion.")
    issorted(t, lt=<=) || error("The time points must be strictly sorted.")
    int64(length(t)) == n || error("Number of time points must be equal to the vector holding the time points.")
    0 < h < 1 || error("Hurst index must be between 0 and 1.")
    new(t, n)
  end
end

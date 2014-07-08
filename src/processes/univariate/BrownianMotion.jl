immutable BrownianMotion <: ContinuousUnivariateStochasticProcess
  timepoints::Vector{Float64}
  npoints::Int64

  function BrownianMotion(t::Vector{Float64}, n::Int64)
    t[1] > 0.0 || error("First provided time point must be positive for Brownian motion.")
    issorted(t, lt=<=) || error("The time points must be strictly sorted.")
    int64(length(t)) == n || error("Number of time points must be equal to the vector holding the time points.")
    new(t, n)
  end
end

BrownianMotion(t::Vector{Float64}) = BrownianMotion(t, int64(length(t)))
BrownianMotion(t::Ranges) = BrownianMotion(collect(t), int64(length(t)))
BrownianMotion(t::Float64, n::Int64) = BrownianMotion(t/n:t/n:t)
BrownianMotion(t::Float64) = BrownianMotion([t], 1)

BrownianMotion(t::Matrix{Float64}) = BrownianMotion[BrownianMotion(t[:, i]) for i = 1:size(t, 2)]
BrownianMotion(t::Ranges, n::Int) = BrownianMotion[BrownianMotion(t) for i = 1:n]
BrownianMotion(t::Float64, npoints::Int64, npaths::Int) = BrownianMotion[BrownianMotion(t, npoints) for i = 1:npaths]

function rand!(p::BrownianMotion, x::Vector{Float64})
  x[1] = rand(Normal(0.0, sqrt(p.timepoints[1])))
  for i = 2:p.npoints
    x[i] = rand(Normal(0.0, sqrt(p.timepoints[i]-p.timepoints[i-1])))+x[i-1]
  end
  x
end

function rand!(p::Vector{BrownianMotion}, x::Matrix{Float64})
  for j = 1:length(p)
    x[1, j] = rand(Normal(0.0, sqrt(p[j].timepoints[1])))
    for i = 2:p[j].npoints
      x[i, j] = rand(Normal(0.0, sqrt(p[j].timepoints[i]-p[j].timepoints[i-1])))+x[i-1, j]
    end
  end
  x
end

rand(p::BrownianMotion) = rand!(p, Array(Float64, p.npoints))

rand(p::Vector{BrownianMotion}) = rand!(p, Array(Float64, p[1].npoints, length(p)))

immutable WienerProcess <: ContinuousUnivariateStochasticProcess
  timepoints::Vector{Float64}
  npoints::Int64

  function WienerProcess(t::Vector{Float64}, n::Int64)
    t[1] > 0.0 || error("Second time point must be positive for a Wiener process.")
    issorted(t, lt=<=) || error("The time points must be strictly sorted.")
    int64(length(t)) == n || error("Number of time points must be equal to the vector holding the time points.")
    new(t, n)
  end
end

WienerProcess(t::Vector{Float64}) = WienerProcess(t, int64(length(t)))
WienerProcess(t::Ranges) = WienerProcess(collect(t), int64(length(t)))
WienerProcess(t::Float64, n::Int64) = WienerProcess(t/n:t/n:t)
WienerProcess(t::Float64) = WienerProcess([t], 1)

WienerProcess(t::Matrix{Float64}) = WienerProcess[WienerProcess(t[:, i]) for i = 1:size(t, 2)]
WienerProcess(t::Ranges, n::Int) = WienerProcess[WienerProcess(t) for i = 1:n]
WienerProcess(t::Float64, npoints::Int64, npaths::Int) = WienerProcess[WienerProcess(t, npoints) for i = 1:npaths]

typealias BrownianMotion WienerProcess

function rand!(p::WienerProcess, x::Vector{Float64})
  x[1] = rand(Normal(0.0, sqrt(p.timepoints[1])))
  for i = 2:p.npoints
    x[i] = rand(Normal(0.0, sqrt(p.timepoints[i]-p.timepoints[i-1])))+x[i-1]
  end
  x
end

function rand!(p::Vector{WienerProcess}, x::Matrix{Float64})
  for j = 1:length(p)
    x[1, j] = rand(Normal(0.0, sqrt(p[j].timepoints[1])))
    for i = 2:p[j].npoints
      x[i, j] = rand(Normal(0.0, sqrt(p[j].timepoints[i]-p[j].timepoints[i-1])))+x[i-1, j]
    end
  end
  x
end

rand(p::WienerProcess) = rand!(p, Array(Float64, p.npoints))

rand(p::Vector{WienerProcess}) = rand!(p, Array(Float64, p[1].npoints, length(p)))

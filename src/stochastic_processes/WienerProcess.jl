immutable WienerProcess <: ContinuousUnivariateStochasticProcess
  timepoints::Vector{Float64}
  npoints::Int

  function WienerProcess(timepoints::Vector{Float64}, npoints::Int)
    timepoints[1] > 0.0 || error("Second time point must be positive for a Wiener process.")
    issorted(timepoints, lt=<=) || error("The time points must be strictly sorted.")
    new(timepoints, length(timepoints))
  end
end

WienerProcess(timepoints::Vector{Float64}) = WienerProcess(timepoints, length(timepoints))
WienerProcess(timepoints::Ranges) = WienerProcess(collect(timepoints), length(timepoints))
WienerProcess(t::Float64, n::Int) = WienerProcess(0:t/n:t, n)
WienerProcess(t::Float64) = WienerProcess([t], 1)

function rand(p::WienerProcess)
  samples = Array(Float64, p.npoints)

  samples[1] = rand(Normal(0.0, sqrt(p.timepoints[1])))
  for i = 2:p.npoints
    samples[i] = rand(Normal(0.0, sqrt(p.timepoints[i]-p.timepoints[i-1])))
  end

  samples
end

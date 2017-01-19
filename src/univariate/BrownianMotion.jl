immutable BrownianMotion <: ContinuousUnivariateStochasticProcess
  t::Vector{Float64}
  n::Int64

  function BrownianMotion(t::Vector{Float64}, n::Int64)
    t[1] == 0. || error("Starting time point must be equal to 0.0.")
    issorted(t, lt=<=) || error("The time points must be strictly sorted.")
    length(t) == n || error("Number of time points must be equal to the vector holding the time points.")
    new(t, n)
  end
end

BrownianMotion(t::Vector{Float64}) = BrownianMotion(t, Int64(length(t)))
BrownianMotion(t::Range) = BrownianMotion(collect(t), Int64(length(t)))
BrownianMotion(t::Float64, n::Int64) = BrownianMotion(0.0:t/n:t-t/n)
BrownianMotion(t::Float64) = BrownianMotion([t], 1)

BrownianMotion(t::Matrix{Float64}) = BrownianMotion[BrownianMotion(t[:, i]) for i = 1:size(t, 2)]
BrownianMotion(t::Range, np::Int) = BrownianMotion[BrownianMotion(t) for i = 1:np]
BrownianMotion(t::Float64, n::Int64, np::Int) = BrownianMotion[BrownianMotion(t, n) for i = 1:np]

function rand!(p::BrownianMotion, x::Vector{Float64})
  x[1] = 0.
  for i = 2:p.n
    x[i] = rand(Normal(0.0, sqrt(p.t[i]-p.t[i-1])))+x[i-1]
  end
  x
end

function rand!(p::Vector{BrownianMotion}, x::Matrix{Float64})
  for j = 1:length(p)
    x[1, j] = 0.
    for i = 2:p[j].n
      x[i, j] = rand(Normal(0.0, sqrt(p[j].t[i]-p[j].t[i-1])))+x[i-1, j]
    end
  end
  x
end

rand(p::BrownianMotion) = rand!(p, Array(Float64, p.n))

rand(p::Vector{BrownianMotion}) = rand!(p, Array(Float64, p[1].n, length(p)))

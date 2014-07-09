immutable FBM <: ContinuousUnivariateStochasticProcess
  timepoints::Vector{Float64}
  npoints::Int64
  hurst::Float64

  function FBM(t::Vector{Float64}, n::Int64, h::Float64)
    t[1] > 0.0 || error("First provided time point must be positive for Brownian motion.")
    issorted(t, lt=<=) || error("The time points must be strictly sorted.")
    int64(length(t)) == n || error("Number of time points must be equal to the vector holding the time points.")
    0 < h < 1 || error("Hurst index must be between 0 and 1.")
    new(t, n, h)
  end
end

FBM(t::Vector{Float64}, h::Float64) = FBM(t, int64(length(t)), h)
FBM(t::Ranges, h::Float64) = FBM(collect(t), int64(length(t)), h)
FBM(t::Float64, n::Int64, h::Float64) = FBM(t/n:t/n:t, h)
FBM(t::Float64, h::Float64) = FBM([t], 1, h)

FBM(t::Matrix{Float64}, h::Float64) = FBM[FBM(t[:, i], h) for i = 1:size(t, 2)]
FBM(t::Ranges, n::Int, h::Float64) = FBM[FBM(t, h) for i = 1:n]
FBM(t::Float64, npoints::Int64, npaths::Int, h::Float64) = FBM[FBM(t, npoints, h) for i = 1:npaths]

function cov(p::FBM, i::Int64, j::Int64)
  twoh::Float64 = 2*p.hurst
  0.5*((p.timepoints[i])^twoh+(p.timepoints[j])^twoh-abs(p.timepoints[i]-p.timepoints[j])^twoh)
end

function cov(p::FBM)
  c = Array(Float64, p.npoints, p.npoints)

  for i = 1:p.npoints
    for j = 1:i
      c[i, j] = cov(p, i, j)
    end
  end

  for i = 1:p.npoints
    for j = (i+1):p.npoints
      c[i, j] = c[j, i]
    end
  end

  c
end

rand_chol!(p::FBM, x::Vector{Float64}) = chol(cov(p), :L)*randn(p.npoints)

rand_chol(p::FBM) = rand!(p, Array(Float64, p.npoints))

function rand_chol(p::Vector{FBM})
  np = length(p)

  if np > 1
    for i = 2:np
      p[1].npoints == p[i].npoints || error("All FBM must have same npoints.")
    end
  end

  x = Array(Float64, p[1].npoints, np)

  for i = 1:np
    x[:, 1] = rand_chol(p[i])
  end

  x
end

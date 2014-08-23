# Fractional Brownian motion
immutable FBM <: ContinuousUnivariateStochasticProcess
  t::Vector{Float64}
  n::Int64
  h::Float64 # Hurst index

  function FBM(t::Vector{Float64}, n::Int64, h::Float64)
    t[1] == 0. || error("Starting time point must be equal to 0.0.")
    issorted(t, lt=<=) || error("The time points must be strictly sorted.")
    int64(length(t)) == n || error("Number of time points must be equal to the vector holding the time points.")
    0. < h < 1. || error("Hurst index must be between 0 and 1.")
    new(t, n, h)
  end
end

FBM(t::Vector{Float64}, h::Float64) = FBM(t, int64(length(t)), h)
FBM(t::Ranges, h::Float64) = FBM(collect(t), int64(length(t)), h)
FBM(t::Float64, n::Int64, h::Float64) = FBM(t/n:t/n:t, h)
FBM(t::Float64, h::Float64) = FBM([t], 1, h)

FBM(t::Matrix{Float64}, h::Float64) = FBM[FBM(t[:, i], h) for i = 1:size(t, 2)]
FBM(t::Ranges, np::Int, h::Float64) = FBM[FBM(t, h) for i = 1:np]
FBM(t::Float64, n::Int64, np::Int, h::Float64) = FBM[FBM(t, n, h) for i = 1:np]

# Fractional Gaussian noise
immutable FGN <: ContinuousUnivariateStochasticProcess
  σ::Float64
  h::Float64 # Hurst index

  function FGN(σ::Float64, h::Float64)
    σ > 0. || error("Standard deviation must be positive.")
    0. < h < 1. || error("Hurst index must be between 0 and 1.")
    new(σ, h)
  end
end

FGN(h::Float64) = FGN(1., h)

convert(::Type{FBM}, p::BrownianMotion) = FBM(p.t, p.n, 0.5)
convert(::Type{FGN}, p::BrownianMotion) = FGN((p.t[end]/(p.n-1))^0.5, 0.5)
convert(::Type{FGN}, p::FBM) = FGN((p.t[end]/(p.n-1))^p.h, p.h)

function autocov(p::FGN, i::Int64, j::Int64)
  twoh::Float64 = 2*p.h
  0.5*abs2(p.σ)*(abs(j-i+1)^twoh+abs(j-i-1)^twoh-2*abs(j-i)^twoh)
end

function autocov!(c::Matrix{Float64}, p::FGN)
  n::Int64 = size(c, 1)

  for i = 1:n
    for j = 1:i
      c[i, j] = autocov(p, i, j)
    end
  end

  for i = 1:n
    for j = (i+1):n
      c[i, j] = c[j, i]
    end
  end

  c
end

autocov(p::FGN, maxlag::Int64) = autocov!(Array(Float64, maxlag, maxlag), p)

function autocov!(y::Vector{Float64}, p::FGN, lags::IntegerVector)
  nlags = length(lags)
  twoh::Float64 = 2*p.h

  for i = 1:nlags
    y[i] = 0.5*abs2(p.σ)*(abs(lags[i]+1)^twoh+abs(lags[i]-1)^twoh-2*abs(lags[i])^twoh)
  end

  y
end

autocov(p::FGN, lags::IntegerVector) = autocov!(Array(Float64, length(lags)), p, lags)

function autocov(p::FBM, i::Int64, j::Int64)
  twoh::Float64 = 2*p.h
  0.5*((p.t[i])^twoh+(p.t[j])^twoh-abs(p.t[i]-p.t[j])^twoh)
end

function autocov!(c::Matrix{Float64}, p::FBM)
  n::Int64 = p.n-1

  for i = 1:n
    for j = 1:i
      c[i, j] = autocov(p, i+1, j+1)
    end
  end

  for i = 1:n
    for j = (i+1):n
      c[i, j] = c[j, i]
    end
  end

  c
end

function autocov(p::FBM)
  n::Int64 = p.n-1
  autocov!(Array(Float64, n, n), p)
end

### rand_chol generates FBM using the method based on Cholesky decomposition.
### T. Dieker, Simulation of Fractional Brownian Motion, master thesis, 2004.
### The complexity of the algorithm is O(n^3), where n is the number of FBM samples.
function rand_chol(p::FBM; fbm::Bool=true)
  w = chol(autocov(p), :L)*randn(p.n-1)

  # If fbm is true return FBM, otherwise return FGN
  w = [0., w]
  if !fbm
    w = diff(w)
  end

  w
end

function rand_chol(p::Vector{FBM}; fbm::Bool=true)
  n::Int64 = fbm ? p[1].n : p[1].n-1
  np = length(p)
  x = Array(Float64, n, np)

  for i = 1:np
    x[:, i] = rand_chol(p[i], fbm=fbm)
  end

  x
end

### rand_fft generates FBM using fast Fourier transform (FFT).
### The time interval of FBM is [0, 1] with a stepsize of 2^p, where p is a natural number.
### The algorithm is known as the Davies-Harte method or the method of circular embedding.
### R.B. Davies and D.S. Harte, Tests for Hurst Effect, Biometrika, 74 (1), 1987, pp. 95-102.
### For a more recent publication, see P.F. Craigmile, Simulating a Class of Stationary Gaussian Processes Using the
### Davies–Harte Algorithm, With Application to Long Memory Processes, Journal of Time Series Analysis, 24 (5), 2003,
### pp. 505-511.
### The complexity of the algorithm is O(n*log(n)), where n=2^p is the number of FBM samples.
function rand_fft(p::FBM; fbm::Bool=true)
  # Determine number of points of simulated FBM
  pnmone::Int64 = p.n-1
  n::Int64 = 2^ceil(log2(pnmone))

  # Compute autocovariance sequence of underlying FGN
  c = Array(Float64, n+1)
  autocov!(c, FGN((1/n)^p.h, p.h), 0:n)

  # Compute square root of eigenvalues of circular autocovariance sequence
  l = real(fft([c, c[end-1:-1:2]]))
  all(i->(i>0), l) || error("Non-positive eigenvalues encountered.")
  lsqrt = sqrt(l)

  # Simulate standard random normal variables
  twon::Int64 = 2*n
  z = randn(twon)

  # Generate fractional Gaussian noise (retain only the first p.n-1 values)
  x = sqrt(0.5)*lsqrt[2:n].*complex(z[2*(2:n)-2], z[2*(2:n)-1])
  w = real(bfft([lsqrt[1]*z[1], x, lsqrt[n+1]*z[twon], conj(reverse(x))]))[1:pnmone]/sqrt(twon)

  # If fbm is true return FBM, otherwise return FGN
  if fbm
    w = [0., cumsum(w)]
  end

  w
end

function rand_fft(p::Vector{FBM}; fbm::Bool=true)
  n::Int64 = fbm ? p[1].n : p[1].n-1
  np = length(p)
  x = Array(Float64, n, np)

  for i = 1:np
    x[:, i] = rand_fft(p[i], fbm=fbm)
  end

  x
end

function rand(p::FBM; fbm::Bool=true, rtype::Symbol=:fft)
  if rtype == :chol
    rand_chol(p; fbm=fbm)
  elseif rtype == :fft
    rand_fft(p; fbm=fbm)
  else
  end
end

function rand(p::Vector{FBM}; fbm::Bool=true, rtype::Symbol=:fft)
  if rtype == :chol
    rand_chol(p; fbm=fbm)
  elseif rtype == :fft
    rand_fft(p; fbm=fbm)
  else
  end
end

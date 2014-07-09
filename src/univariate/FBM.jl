immutable FBM <: ContinuousUnivariateStochasticProcess
  timepoints::Vector{Float64}
  npoints::Int64
  hurst::Float64

  function FBM(t::Vector{Float64}, n::Int64, h::Float64)
    t[1] == 0.0 || error("Starting time point must be equal to 0.0.")
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
  npoints::Int64 = p.npoints-1
  c = Array(Float64, npoints, npoints)

  for i = 1:npoints
    for j = 1:i
      c[i, j] = cov(p, i+1, j+1)
    end
  end

  for i = 1:npoints
    for j = (i+1):npoints
      c[i, j] = c[j, i]
    end
  end

  c
end

### rand_chol generates FBM using the method based on Cholesky decomposition.
### T. Dieker, Simulation of Fractional Brownian Motion, master thesis, 2004.
### The complexity of the algorithm is O(n^3), where n is the number of FBM samples.
rand_chol(p::FBM) = [0., chol(cov(p), :L)*randn(p.npoints-1)]

function rand_chol(p::Vector{FBM})
  np::Int64 = length(p)

  if np > 1
    for i = 2:np
      p[1].npoints == p[i].npoints || error("All FBM must have same number of points.")
    end
  end

  x = Array(Float64, p[1].npoints, np)

  for i = 1:np
    x[:, i] = rand_chol(p[i])
  end

  x
end

### rand_fft generates FBM using fast Fourier transform (FFT).
### The time interval of FBM is [0, 1] with a stepsize of 2^p, where p is a natural number.
### The algorithm is known as the Davis-Harte method or the method of circular embedding.
### R.B. Davies and D.S. Harte, Tests for Hurst Effect, Biometrika, 74 (1987), pp. 95–102.
### The complexity of the algorithm is O(n*log(n)), where n=2^p is the number of FBM samples.
function rand_fft(p::FBM)
  # Determine number of points of simulated FBM
  log2n = log2(p.npoints-1)
  npoints::Int64 = isinteger(log2n) ? p.npoints : 2^floor(log2n)

  # Construct circular covariant matrix
  c = Array(Float64, npoints+1)
  twop::Float64 = 2*p.hurst
  c[1] = 1
  for k = 1:npoints
    c[k+1] = 0.5*((k+1)^twop+(k-1)^twop-2*k^twop)
  end
  c = [c, c[end-1:-1:2]]

  # Compute the eigenvalues of the circular covariant matrix
  twonpoints = 2*npoints
  l = real(fft(c))/(twonpoints)

  # Derive Fractional Gaussian Noise (FGN)
  w = fft(sqrt(l).*complex(randn(twonpoints), randn(twonpoints)))
  w = npoints^(-p.hurst)*real(W[1:twonpoints+1])
end

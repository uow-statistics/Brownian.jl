# Fractional Brownian motion
struct FBM <: ContinuousUnivariateStochasticProcess
    t::Vector{Float64}
    n::Int64
    h::Float64 # Hurst index

    function FBM(t::Vector{Float64}, n::Int64, h::Float64)
        t[1] == 0.0 || error("Starting time point must be equal to 0.0.")
        issorted(t, lt = <=) || error("The time points must be strictly sorted.")
        length(t) == n || error(
            "Number of time points must be equal to the vector holding the time points.",
        )
        0 < h < 1 || error("Hurst index must be between 0 and 1.")
        new(t, n, h)
    end
end

FBM(t::Vector{Float64}, h::Float64) = FBM(t, Int64(length(t)), h)
FBM(t::AR, h::Float64) where {AR<:AbstractRange} = FBM(collect(t), Int64(length(t)), h)
FBM(t::Float64, n::Int64, h::Float64) = FBM(0.0:t/n:t-t/n, h)
FBM(t::Float64, h::Float64) = FBM([t], 1, h)

FBM(t::Matrix{Float64}, h::Float64) = FBM[FBM(t[:, i], h) for i = 1:size(t, 2)]
FBM(t::AR, np::Int, h::Float64) where {AR<:AbstractRange} = FBM[FBM(t, h) for i = 1:np]
FBM(t::Float64, n::Int64, np::Int, h::Float64) = FBM[FBM(t, n, h) for i = 1:np]

# Fractional Gaussian noise
struct FGN <: ContinuousUnivariateStochasticProcess
    σ::Float64
    h::Float64 # Hurst index

    function FGN(σ::Float64, h::Float64)
        σ > 0.0 || error("Standard deviation must be positive.")
        0.0 < h < 1.0 || error("Hurst index must be between 0 and 1.")
        new(σ, h)
    end
end

FGN(h::Float64) = FGN(1.0, h)

convert(::Type{FBM}, p::BrownianMotion) = FBM(p.t, p.n, 0.5)
convert(::Type{FGN}, p::BrownianMotion) = FGN((p.t[end] / (p.n - 1))^0.5, 0.5)
convert(::Type{FGN}, p::FBM) = FGN((p.t[end] / (p.n - 1))^p.h, p.h)

function autocov(p::FGN, i::Int64, j::Int64)
    twoh::Float64 = 2 * p.h
    0.5 * abs2(p.σ) * (abs(j - i + 1)^twoh + abs(j - i - 1)^twoh - 2 * abs(j - i)^twoh)
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

autocov(p::FGN, maxlag::Int64) = autocov!(Array{Float64}(undef, maxlag, maxlag), p)

function autocov!(y::Vector{Float64}, p::FGN, lags::AbstractVector{<:Integer})
    nlags = length(lags)
    twoh::Float64 = 2 * p.h

    for i = 1:nlags
        y[i] =
            0.5 *
            abs2(p.σ) *
            (abs(lags[i] + 1)^twoh + abs(lags[i] - 1)^twoh - 2 * abs(lags[i])^twoh)
    end

    y
end

autocov(p::FGN, lags::AbstractVector{<:Integer}) =
    autocov!(Array{Float64}(undef, length(lags)), p, lags)

function autocov(p::FBM, i::Int64, j::Int64)
    twoh::Float64 = 2 * p.h
    0.5 * ((p.t[i])^twoh + (p.t[j])^twoh - abs(p.t[i] - p.t[j])^twoh)
end

function autocov!(c::Matrix{Float64}, p::FBM)
    n::Int64 = p.n - 1

    for i = 1:n
        for j = 1:i
            c[i, j] = autocov(p, i + 1, j + 1)
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
    n::Int64 = p.n - 1
    autocov!(Array{Float64}(undef, n, n), p)
end

### rand_chol generates FBM using the method based on Cholesky decomposition.
### T. Dieker, Simulation of Fractional Brownian Motion, master thesis, 2004.
### The complexity of the algorithm is O(n^3), where n is the number of FBM samples.
function rand_chol(p::FBM, fbm::Bool = true)

    w = (cholesky(autocov(p)).U') * randn(p.n - 1)

    # If fbm is true return FBM, otherwise return FGN
    insert!(w, 1, 0.0)
    if !fbm
        w = diff(w)
    end

    w
end

### rand_fft generates FBM using fast Fourier transform (FFT).
### The time interval of FBM is [0, 1] with a stepsize of 2^p, where p is a natural number.
### The algorithm is known as the Davies-Harte method or the method of circular embedding.
### R.B. Davies and D.S. Harte, Tests for Hurst Effect, Biometrika, 74 (1), 1987, pp. 95-102.
### For a more recent publication, see P.F. Craigmile, Simulating a Class of Stationary Gaussian Processes Using the
### Davies–Harte Algorithm, With Application to Long Memory Processes, Journal of Time Series Analysis, 24 (5), 2003,
### pp. 505-511.
### The complexity of the algorithm is O(n*log(n)), where n=2^p is the number of FBM samples.
function rand_fft(p::FBM, fbm::Bool = true)
    # Determine number of points of simulated FBM
    pnmone = p.n - 1
    n = 1 << Int(ceil(log2(pnmone)))

    # Compute autocovariance sequence of underlying FGN
    c = Vector{Float64}(undef, n + 1) # old was Array{Float64}(undef, n+1)
    autocov!(c, FGN((1 / n)^p.h, p.h), 0:n)

    # Compute square root of eigenvalues of circular autocovariance sequence
    l = real(fft(cat(c, c[end-1:-1:2], dims = (1))))
    all(i -> (i > 0), l) || error("Non-positive eigenvalues encountered.")
    lsqrt = sqrt.(l)

    # Simulate standard random normal variables
    twon::Int64 = 2 * n
    z = randn(twon)

    # Generate fractional Gaussian noise (retain only the first p.n-1 values)
    x =
        sqrt(0.5) * lsqrt[2:n] .*
        complex.(z[collect(2 .* (2:n) .- 2)], z[collect(2 .* (2:n) .- 1)]) #2*(2:n)-2
    w =
        real(
            bfft(
                cat(lsqrt[1] * z[1], x, lsqrt[n+1] * z[twon], conj(reverse(x)), dims = (1)),
            ),
        )[1:pnmone] / sqrt(twon)

    # If fbm is true return FBM, otherwise return FGN
    if fbm
        w = insert!(cumsum(w), 1, 0.0)
    end

    w
end

### Implementation of "exact discrete" method for simulating Riemann-Liouville fBm
### Based on Muniandy, S. & Lim, S. Modeling of locally self-similar processes using multifractional
### Brownian motion of Riemann-Liouville type. Physical Review E 63, 046104. ISSN: 1063-651X (2001).
### Specifically, Eqn. 17, 18 and 19.
### From Muniandy (2001), Reimann-Liouville fBm is defined as,
### B_H(t)=\frac{1}{\Gamma(H+0.5)}\int^{t}_{0}(t-s)^(H-0.5)dB(s), t\geq0.
### Here we focus on the discrete time t_j=j\Delta t approximation,
### B_H(t_j)=\frac{1}{\Gamma(H+0.5)}\sum^{j}_{i=1}\int^{i\Delta t}_{(i-1)\Delta t}(t_j-\tau)^(H-0.5)dB(\tau)
### The exact solution to the interior integral results in a weighting function, implemented here with the :Exact key.
### An "improved" weighting function was proposed by
### Rambaldi, S. & Pinazza, O. An accurate fractional Brownian motion generator. Physica A 208, 21–30 (1994).
### We have implemented this here with the :Improved key.
function rand_rl(p::FBM, fbm::Bool = true; wtype::Symbol = :exact)
    # We are going to require that the time points are evenly spaced to keep things simple
    dt = unique(diff(p.t))
    if length(dt) != 1
        error(
            "For this simulation technique, the fBm process must be sampled across a uniform temporal grid",
        )
    end

    dt = dt[1]

    w_mat = zeros(p.n - 1, p.n - 1) # Prealocate a weight matrix. Note that this will be upper triangular

    for j = 1:p.n-1
        for i = 1:j
            w_mat[i, j] = weight(p.h, (j - i + 1) * dt, dt, wtype)
        end
    end

    # Multiply a vector of white noise by the weight matrix and scale appropriately.
    X = dropdims(sqrt(2) * (randn(p.n - 1)' * w_mat), dims = (1)) .* sqrt(dt)

    insert!(X, 1, 0.0)  # Snap to zero. Not clear if this causes a discontinuity in the correlation structure.
    #It is worth considering alternative for the above line: generate path X of length p.n, then let X = X-X[1].

    return X
end

weight(H::Float64, t::Float64, dt::Float64, method::Symbol = :exact) =
    weight(H, t, dt, Val{method})

function weight(H::Float64, t::Float64, dt::Float64, ::Type{Val{:exact}})
    nom = t^(H + 0.5) - (t - dt)^(H + 0.5)
    denom = gamma(H + 0.5) * (H + 0.5) * dt
    return nom / denom
end

function weight(H::Float64, t::Float64, dt::Float64, ::Type{Val{:improved}})
    nom = t^(2 * H) - (t - dt)^(2 * H)
    denom = 2 * H * dt
    return sqrt(nom / denom) / gamma(H + 0.5)
end

### API functions

function rand(p::FBM; fbm::Bool = true, method::Symbol = :fft, args...)
    if method == :fft
        rand_fft(p, fbm)
    elseif method == :rl
        rand_rl(p, fbm; args...)
    elseif method == :chol
        rand_chol(p, fbm)
    else
        error("Non-existing method $method")
    end
end

function rand(p::Vector{FBM}; fbm::Bool = true, method::Symbol = :fft, args...)
    n::Int64 = fbm ? p[1].n : p[1].n - 1
    np = length(p)
    x = Array{Float64}(undef, n, np)

    for i = 1:np
        x[:, i] = rand(p[i]; fbm = fbm, method = method, args...)
    end

    x
end

function chol_update(r0::Array{Float64,2}, A::Array{Float64,2})
    # r0 is a lower triangular matrix found via cholesky decomposition of
    # A(1:end-1,1:end-1), therefore r0 is a matrix of size (N-1)-by-(N-1) when A
    # is of size N-by-N. A is a positive symetric semidefinite matrix by construction.

    m = size(A, 1) # Assumes A is square
    r1 = zeros(m, m)
    r1[1:end-1, 1:end-1] = r0

    @inbounds for i = 1:m-1
        r1[m, i] = (1.0 / r1[i, i]) * (A[m, i]-sum(r1[m, 1:i] .* r1[i, 1:i]))[1]
    end

    r1[m, m] = sqrt(A[m, m] - sum(r1[m, 1:m-1] .* r1[m, 1:m-1]))[1]

    return r1
end

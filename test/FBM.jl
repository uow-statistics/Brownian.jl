@test convert(FGN, FBM(0:0.5:10, 0.4)) == FGN(0.757858283255199, 0.4)

p = FBM(0:1/2^10:1, 0.4)

# fBm using FFT approach

rand(p)
rand(p, fbm = false)
rand([p, p])

# fBm using Riemann-Liouville approach

rand(p, method = :rl)
rand(p, method = :rl, wtype = :improved)
rand([p, p], method = :rl)
rand([p, p], method = :rl, wtype = :improved)

# fBm using Cholesky approach

rand(p, method = :chol)
rand(p, fbm = false, method = :chol)
rand([p, p], method = :chol)

# Checking function for efficient Cholesky update

p = FBM(0:0.1:0.1*100, 0.4)
q = FBM(0:0.1:0.1*101, 0.4)

P = autocov(p)
Q = autocov(q)

using LinearAlgebra
c = cholesky(P).U'
chol_update(convert(Array{Float64,2}, c), Q)

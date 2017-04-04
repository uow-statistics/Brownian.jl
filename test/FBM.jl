using Brownian
using Base.Test

@test convert(FGN, FBM(0:0.5:10, 0.4)) == FGN(0.757858283255199,0.4)

p = FBM(0:1/2^10:1, 0.4)

# fBm using FFT approach

rand(p)
rand(p, fbm=false)
rand([p, p])

# fBm using Riemann-Liouville approach

rand(p, method=:rl)
rand(p, method=:rl, wtype=:improved)

# fBm using Cholesky approach

rand(p, method=:chol)
rand(p, fbm=false, method=:chol)
rand([p, p], method=:chol)

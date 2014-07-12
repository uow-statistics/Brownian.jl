using Brownian
using Base.Test

@test convert(FGN, FBM(0:0.5:10, 0.4)) == FGN(0.757858283255199,0.4)

p = FBM(0:1/2^10:1, 0.4)

rand(p)
rand(p, fbm=false)
rand([p, p])
rand(p, rtype=:chol)
rand(p, fbm=false, rtype=:chol)
rand([p, p], rtype=:chol)

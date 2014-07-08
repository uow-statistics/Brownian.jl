immutable FGN <: ContinuousUnivariateStochasticProcess
  hurst::Float64

  function FGN(h::Float64)
    0 < h < 1 || error("Hurst index must be between 0 and 1.")
    new(h)
  end
end

function cov(p::FGN, i::Int64, j::Int64)
  twoh::Float64 = 2*p.hurst
  0.5*(abs(j-i+1)^twoh+abs(j-i-1)^twoh-2*abs(j-i)^twoh)
end

function cov(p::FGN, n::Int64)
  c = Array(Float64, n, n)

  for i = 1:n
    for j = 1:i
      c[i, j] = cov(p, i, j)
    end
  end

  for i = 1:n
    for j = (i+1):n
      c[i, j] = c[j, i]
    end
  end

  c
end

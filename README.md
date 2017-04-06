## Brownian

[![Build Status](https://travis-ci.org/UniversityofWarwick/Brownian.jl.png)](https://travis-ci.org/UniversityofWarwick/Brownian.jl)

The Julia `Brownian` package is aimed at providing a centralized repository of algorithms for simulating Brownian-based
stochastic processes. More precisely, the package currently provides routines for random sampling from
* one-dimensional Brownian motion via random walk,
* one-dimensional fractional Brownian motion (FBM) and one-dimensional fractional Gaussian noise (FGN) via the Cholesky
decomposition method or the Davies-Harte method, which makes use of fast Fourier transforms,
* one-dimensional Riemann-Liouville fractional Brownian motion (FBM) via an exact discrete method.

The future roadmap would be to provide implementations for sampling from
* one-dimensional Brownian motion via Brownian bridge and via multivariate normals,
* one-dimensional fractional Brownian motion using the Hosking method,
* multidimensional Brownian and fractional Brownian motion,
* reflected Brownian motion (RBM).

Willing developers are welcome to contribute to the package.

### Short tutorial

### Example 1: simulation of Brownian motion

To simulate Brownian motion at the time points (0, 0.1, 0.5, 0.75, 1), use the following snippet:

```
using Brownian

p = BrownianMotion([0, 0.1, 0.5, 0.75, 1])

rand(p)
```

### Example 2: simulation of FBM and FGN

Suppose that interest is in simulating fractional Brownian motion with Hurst index equal to 0.4 in the time interval
[0, 1] with a time step of 1/2^n for some natural n (for example n=10).

```
using Brownian

p = FBM(0:1/2^10:1, 0.4)

# Using the Davies-Harte algorithm, which relies on fast Fourier transforms (FFT)
rand(p)

# Using a method based on the Cholesky decomposition of the covariance matrix of FBM
rand(p, method=:chol)

# Using an exact discrete method for simulating Riemann-Liouville FBM
rand(p, method=:rl)
```

To simulate fractional Gaussian noise with the same Hurst index,

```
# Using the Davies-Harte algorithm
rand(p, fbm=false)

# Using the Cholesky method
rand(p, fbm=false, method=:chol)
```

Note that fractional Brownian motion is obtained from fractional Gaussian noise by taking cumulative sums (and
conversely FGN is computed from FBM by differencing).

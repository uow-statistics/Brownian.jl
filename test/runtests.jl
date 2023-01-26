using Test, Brownian
println("Running tests:")

tests = ["BrownianMotion", "FBM"]

for t in tests
    test_fn = "$t.jl"
    println("  * $test_fn")
    include(test_fn)
end

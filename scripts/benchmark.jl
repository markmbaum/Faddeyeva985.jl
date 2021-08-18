using Faddeyeva985, SpecialFunctions, BenchmarkTools

##

#computes "exact" faddeyeva
wofz(z::Complex)::Complex = erfcx(-im*z)

grids = [
    (-500, 500, 4000,  1e-5, 1e5, 1000),
    (-200, 200, 4000, 1e-20, 1e4, 1000),
    (-10,   10, 4000, 1e-20, 1e4, 1000),
    (-500, 500, 1000, 1e-20, 1e5, 2500)
]

meantime(b) = mean(b.times)/1e9

##

for (n,grid) ∈ enumerate(grids)
    x₁, x₂, nx, y₁, y₂, ny = grid
    x = LinRange(x₁, x₂, nx)
    y = 10 .^ LinRange(log10(y₁), log10(y₂), ny)
    X = x' .* ones(length(y))
    Y = y .* ones(length(x))'
    Z = X .+ im*Y
    println("\ngrid $n")
    #full complex evaluation
    b₁ = meantime(@benchmark wofz.($Z))
    println("  exact complex result: $b₁ seconds")
    b₂ = meantime(@benchmark faddeyeva.($Z))
    println("  apprx complex result: $b₂ seconds")
    println("    $(b₁/b₂)")
    #real part only
    b₁ = meantime(@benchmark real.(wofz.($Z)))
    println("  exact real result: $b₁ seconds")
    b₂ = meantime(@benchmark faddeyeva.($X, $Y))
    println("  apprx real result: $b₂ seconds")
    println("    $(b₁/b₂)")
end
# Faddeyeva985.jl

[![Build Status](https://github.com/markmbaum/Faddeyeva985.jl/workflows/CI/badge.svg)](https://github.com/markmbaum/Faddeyeva985.jl/actions)
[![Coverage](https://codecov.io/gh/markmbaum/Faddeyeva985.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/markmbaum/Faddeyeva985.jl)

This package exports a single function approximating the [Faddeyeva Function](https://en.wikipedia.org/wiki/Faddeeva_function) (sometimes called Faddeeva). [The algorithm](https://dl.acm.org/doi/10.1145/3119904) was developed by [Mofreh R. Zaghloul](https://cos.uaeu.ac.ae/en/departments/physics/profile.shtml?email=m.zaghloul) at the United Arab Emirates University and is translated to Julia here. It achieves relative accuracy better than 4e-5 over a wide range of arguments for the real and imaginary parts of the result.

The function has two methods. One takes a complex argument and returns the full complex result.
```julia
julia> faddeyeva(1 + 2im)
0.2184926141618737 + 0.09299780803270477im
```
The other takes a complex argument split into its parts and returns only the real part of the result. This spares some calculations.
```julia
julia> faddeyeva(1, 2)
0.2184926141618737
```

A highly accurate version of the function is available through [`erfcx`](https://specialfunctions.juliamath.org/stable/functions_list/#SpecialFunctions.erfcx) in [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), but it's slower.
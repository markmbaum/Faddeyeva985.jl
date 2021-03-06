# Faddeyeva985.jl

[![Build Status](https://github.com/markmbaum/Faddeyeva985.jl/workflows/CI/badge.svg)](https://github.com/markmbaum/Faddeyeva985.jl/actions)
[![codecov](https://codecov.io/gh/markmbaum/Faddeyeva985.jl/branch/main/graph/badge.svg?token=uOEvwf0hpm)](https://codecov.io/gh/markmbaum/Faddeyeva985.jl)

This package exports a single function approximating the [Faddeyeva Function](https://en.wikipedia.org/wiki/Faddeeva_function) (sometimes called Faddeeva). [The algorithm](https://dl.acm.org/doi/10.1145/3119904) was developed by [Mofreh R. Zaghloul](https://cos.uaeu.ac.ae/en/departments/physics/profile.shtml?email=m.zaghloul) at the United Arab Emirates University and is translated to Julia here.

The published description states a relative accuracy better than 4e-5 over wide ranges of inputs. My testing indicates relative error is *slightly* greater than 4e-5 at some points in the same ranges, but never greater than 5e-5.

The Faddeyeva function arises in a handful of physical problems. For example, it is needed to evaluate the [Voigt profile](https://en.wikipedia.org/wiki/Voigt_profile) for molecular absorption lines in radiative transfer applications/modeling.

If you use it for research, cite the original paper:

* Mofreh R. Zaghloul. 2017. Algorithm 985: Simple, Efficient, and Relatively Accurate Approximation for the Evaluation of the Faddeyeva Function. *ACM Trans. Math. Softw.* 44, 2, Article 22 (October 2017), 9 pages. DOI: [https://doi.org/10.1145/3119904](https://doi.org/10.1145/3119904)

-----

The module can be installed from Julia's general registry with `]add Faddeyeva985`

The function has two methods. One takes a complex argument and returns the full complex result.
```julia
julia> using Faddeyeva985
julia> faddeyeva(1 + 2im)
0.2184926141618737 + 0.09299780803270477im
```
The other takes a complex argument split into its parts and returns only the real part of the result. This spares some calculations.
```julia
julia> faddeyeva(1, 2)
0.2184926141618737
```

A highly accurate version of the function is available through [`erfcx`](https://specialfunctions.juliamath.org/stable/functions_list/#SpecialFunctions.erfcx) in [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), but it's slower.

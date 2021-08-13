using Faddeyeva985
using SpecialFunctions
using Test

#computes "exact" real part of faddeyeva
wofz(z::Complex)::Complex = erfcx(-im*z)

#compute the maximum relative errors over a grid
function maxrelerr(x₁, x₂, nx, y₁, y₂, ny)
    mie = 0.0 #maximum relative error of imaginary part
    mre = 0.0 #maximum relative error of real part
    #grid points
    xg = LinRange(x₁, x₂, nx)
    yg = 10 .^ LinRange(log10(y₁), log10(y₂), ny)
    for x ∈ xg, y ∈ yg
        #complex argument
        z = x + im*y
        #exact complex result
        w = wofz(z)
        #985 complex result
        f = faddeyeva(z)
        #relative error of imaginary part
        ie = abs(imag(w) - imag(f))/abs(imag(w))
        #relative error of real part
        re = abs(real(w) - real(f))/abs(real(w))
        #store the highest values
        if ie > mie
            mie = ie
        end
        if re > mre
            mre = re
        end
    end
    return mie, mre
end

#these test grids are based on those in Table 2 of the original publication: https://doi.org/10.1145/3119904
grids = [
    (-500, 500, 40000,  1e-5, 1e5, 10000),
    (-200, 200, 40000, 1e-20, 1e4, 10000),
    (-10,   10, 40000, 1e-20, 1e4, 10000),
    (-500, 500, 10000, 1e-20, 1e5, 25000)
]

@testset "Faddeyeva985" begin
    for (n,grid) ∈ enumerate(grids)
        imagerr, realerr = maxrelerr(grid...)
        println("grid $n\n  max imag error = $(imagerr*100) %\n  max real error = $(realerr*100) %")
        @test imagerr < 1e-4
        @test realerr < 1e-4
    end
end
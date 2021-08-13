module Faddeyeva985

export faddeyeva

const Θ = 1/√π

const α₀ = 122.60793
const α₁ = 214.38239
const α₂ = 181.92853
const α₃ = 93.15558
const α₄ = 30.180142
const α₅ = 5.9126262
const α₆ = 1/√π

const β₀ = 122.60793
const β₁ = 352.73063
const β₂ = 457.33448
const β₃ = 348.70392
const β₄ = 170.35400
const β₅ = 53.992907
const β₆ = 10.479857

const γ₀ = 36183.31
const γ₁ = 3321.99
const γ₂ = 1540.787
const γ₃ = 219.031
const γ₄ = 35.7668
const γ₅ = 1.320522
const γ₆ = 1/√π

const λ₀ = 32066.6
const λ₁ = 24322.84
const λ₂ = 9022.228
const λ₃ = 2186.181
const λ₄ = 364.2191
const λ₅ = 61.57037
const λ₆ = 1.841439

const s₀ = 38000.0
const s₁ = 256.0
const s₂ = 62.0
const s₃ = 30.0
const t₃ = 1.0e-13
const s₄ = 2.5
const t₄ = 5.0e-9
const t₅ = 0.072

#region 4: Laplace continued fractions, 4 convergents
function region4(z::Complex, x, y)::Complex
    z² = z^2
    (Θ*(-y + im*x))*(z² - 2.5)/(z²*(z² - 3.0) + 0.75)
end

#region 5: Humlicek's w4 (Region 4), part a
function region5a(z::Complex, x²)::Complex
    z² = z^2
    r = γ₀ + z²*(γ₁ + z²*(γ₂ + z²*(γ₃ + z²*(γ₄ + z²*(γ₅ + z²*γ₆)))))
    t = λ₀ + z²*(λ₁ + z²*(λ₂ + z²*(λ₃ + z²*(λ₄ + z²*(λ₅ + z²*(λ₆ + z²))))))
    exp(-x²) + (im*z*r/t)
end

#region 5: Humlicek's w4 (Region 4), part b
function region5b(z::Complex)::Complex
    z² = z^2
    r = γ₀ + z²*(γ₁ + z²*(γ₂ + z²*(γ₃ + z²*(γ₄ + z²*(γ₅ + z²*γ₆)))))
    t = λ₀ + z²*(λ₁ + z²*(λ₂ + z²*(λ₃ + z²*(λ₄ + z²*(λ₅ + z²*(λ₆ + z²))))))
    exp(-z²) + (im*z*r/t)
end

#region 6: Hui's p-6 Approximation
function region6(x, y)::Complex
    q = y - im*x
    r = α₀ + q*(α₁ + q*(α₂ + q*(α₃ + q*(α₄ + q*(α₅ + q*α₆)))))
    t = β₀ + q*(β₁ + q*(β₂ + q*(β₃ + q*(β₄ + q*(β₅ + q*(β₆ + q))))))
    r/t
end

function faddeyeva(z::Complex)::Complex

    x = real(z)
    y = imag(z)
    x² = x*x
    y² = y*y
    s = x² + y²

    #region 1: Laplace continued fractions, 1 convergent
    if s >= s₀
        return (y + im*x)*Θ/s
    end

    #region 2: Laplace continued fractions, 2 convergents
    if s >= s₁
        a = y*(0.5 + s)
        b = x*(s - 0.5)
        d = s^2 + (y² - x²) + 0.25
        return (a + im*b)*(Θ/d)
    end

    #region 3: Laplace continued fractions, 3 convergents
    if s >= s₂
        q = y² - x² + 1.5
        r = 4.0*x²*y²
        a = y*((q - 0.5)*q + r + x²)
        b = x*((q - 0.5)*q + r - y²)
        d = s*(q*q + r)
        return Θ*(a + im*b)/d
    end

    #region 4: Laplace continued fractions, 4 convergents
    if s >= s₃ && y² >= t₃
        return region4(z, x, y)
    end

    #region 5: Humlicek's w4 (Region 4)
    if s > s₄ && y² < t₄
        return region5a(z, x²)
    elseif s > s₄ && y² < t₅
        return region5b(z)
    end

    #region 6: Hui's p-6 Approximation
    return region6(x, y)

end

#real arguments represent z = x + im*y and return only the real part
function faddeyeva(x, y)

    x² = x*x
    y² = y*y
    s = x² + y²

    #------------------------------------------------------
    #regions 1-3 are optimized/simplified for real-only return

    #region 1: Laplace continued fractions, 1 convergent
    if s >= s₀ 
        return y*Θ/s
    end

    #region 2: Laplace continued fractions, 2 convergents
    if s >= s₁
        return y*(0.5 + s)*(Θ/((s^2 + (y² - x²)) + 0.25))
    end

    #region 3: Laplace continued fractions, 3 convergents
    if s >= s₂
        q = y² - x² + 1.5
        r = 4.0*x²*y²
        return Θ*(y*((q - 0.5)*q + r + x²))/(s*(q*q + r))
    end

    #----------------------------------------------------------
    #regions 4+ are not optimized for real-only return

    if s >= s₃ && y² >= t₃
        return real(region4(x + im*y, x, y))
    end

    if s > s₄ && y² < t₄
        return real(region5a(x + im*y, x²))
    elseif s > s₄ && y² < t₅
        return real(region5b(x + im*y))
    end

    return real(region6(x, y))

end

end

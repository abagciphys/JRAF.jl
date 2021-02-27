###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                           SPECİAL FUNCTIONS                                             #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
export LegendreP, NLegendreP, LegendrePG, LegendreQ, NLegendreQ, LaguerreL, LaguerreG, SphericalHarmonicsY,
SphericalHarmonicsYG, SphericalHarmonicsS, SphericalHarmonicsSG, SphericalSpinors, RadoslawSSpinors, RotaD, 
Rotad, RotaT, SlaterORadial, SlaterOrbitalY, SlaterOrbitalS, ChargeDOneC, OneCenterP, SlaterSRadial, 
Hypergeometric0F0, Hypergeometric1F0, Hypergeometric0F1, Hypergeometric1F1, HypergeometricU, WhittakerM, 
WhittakerW, Hypergeometric2F1, AuxiliaryGk, AuxiliaryGm, AuxiliaryGl, AuxiliaryGh, AuxiliaryG, 
AuxiliaryGrlm, AuxiliaryGrlp, AuxiliaryGrh
###########################################################################################################
########################################## LEGENDRE POLYNOMIALS ###########################################(43)
################# Condon-Shortley phases.
################# First kind associated Legendre polynomials.
################# Type 0 and type 1 respectively correspond to type 2 and type 3 in Mathematica and mpmath.
################# Mathematica formoat: LegendreP[n,m,type,x]
function LegendreP(n::acb, m::acb, x::acb, t::UInt)
    z = parent(x)()
    ccall((:acb_hypgeom_legendre_p, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, UInt, Int),
    z, n, m, x, t, parent(n).prec)
    return z
end
function LegendreP(n::acb, m::acb, x::acb, t::Int)
    LegendreP(n, m, x, UInt(t))
end
function LegendreP(n::Complex, m::Complex, x::Complex, t::Int)
    LegendreP(CF(real(n), imag(n)), CF(real(m), imag(m)), CF(real(x), imag(x)), UInt(t))
end
function LegendreP(n::arb, m::arb, x::arb, t::Int)
    LegendreP(CF(n), CF(m), CF(x), t)
end
function LegendreP(n::Real, m::Real, x::Real, t::Int)
    LegendreP(CF(n), CF(m), CF(x), t)
end
###########################################################################################################
################# First kind normalized associated Legendre polynomials
function NLegendreP(n::acb, m::acb, x::acb, t::UInt)
    Sqrt(((2n+1)//2) * (Factorial(n-m)//Factorial(n+m)))*LegendreP(n, m, x, t)
end
function NLegendreP(n::acb, m::acb, x::acb, t::Int)
    NLegendreP(n, m, x, UInt(t))
end
function NLegendreP(n::Complex, m::Complex, x::Complex, t::Int)
    NLegendreP(CF(real(n), imag(n)), CF(real(m), imag(m)), CF(real(x), imag(x)), UInt(t))
end
function NLegendreP(n::arb, m::arb, x::arb, t::Int)
    NLegendreP(CF(n), CF(m), CF(x), t)
end
function NLegendreP(n::Real, m::Real, x::Real, t::Int)
    NLegendreP(CF(n), CF(m), CF(x), t)
end
###########################################################################################################
################# LegendrePG is a Guseinov's formulations for normalized associated Legendre polynomials
################# ((-1)^m) for Condon-Shortley phase
function LegendrePG(n::Int, m::Int, x::acb)
    lowl = Zero(0)
    upl = Int(floor((n-abs(m))/2))

    res1 = Sqrt((2*n+1)/2)*((1-x^2)^(Abs(m)/2))
    res2 = zeros(CF, upl - lowl + 1)

    for k in lowl : upl
        res2[k - lowl + 1] = SphCB(k, n, abs(m))*(x^(n-abs(m)-2*k))
    end
    cs = m < 0 ? RF(1) : Power(-1,m) #For Condon-Shortley Phases#
    return cs * res1 * sum(res2)
end
###########################################################################################################
################# Second kind associated Legendre polynomials
function LegendreQ(n::acb, m::acb, x::acb, t::UInt)
    z = parent(x)()
    ccall((:acb_hypgeom_legendre_q, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, UInt, Int),
    z, n, m, x, t, parent(n).prec)
    return z
end

function LegendreQ(n::acb, m::acb, x::acb, t::Int)
    LegendreQ(n, m, x, UInt(t))
end
function LegendreQ(n::Complex, m::Complex, x::Complex, t::Int)
    LegendreQ(CF(real(n), imag(n)), CF(real(m), imag(m)), CF(real(x), imag(x)), UInt(t))
end
function LegendreQ(n::arb, m::arb, x::arb, t::Int)
    LegendreQ(CF(n), CF(m), CF(x), t)
end
function LegendreQ(n::Real, m::Real, x::Real, t::Int)
    LegendreQ(CF(n), CF(m), CF(x), t)
end
###########################################################################################################
################# Second kind normalized associated Legendre polynomials
function NLegendreQ(n::acb, m::acb, x::acb, t::UInt)
    Sqrt(((2n+1)//2) * (Factorial(n-m)//Factorial(n+m)))*LegendreP(n, m, x, t)
end
function NLegendreQ(n::acb, m::acb, x::acb, t::Int)
    NLegendreQ(n, m, x, UInt(t))
end
function NLegendreQ(n::Complex, m::Complex, x::Complex, t::Int)
    NLegendreQ(CF(real(n), imag(n)), CF(real(m), imag(m)), CF(real(x), imag(x)), UInt(t))
end
function NLegendreQ(n::arb, m::arb, x::arb, t::Int)
    NLegendreQ(CF(n), CF(m), CF(x), t)
end
function NLegendreQ(n::Real, m::Real, x::Real, t::Int)
    NLegendreQ(CF(n), CF(m), CF(x), t)
end
###########################################################################################################
########################################## LAGUERRE POLYNOMIALS ###########################################(43)
################# Associated Laguerre polynomials L_{p}^{q}(x)
function LaguerreL(q::acb, p::acb, x::acb)
    z = parent(x)()
    ccall((:acb_hypgeom_laguerre_l, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Int),
    z, q, p, x, parent(x).prec)
    return z
end
function LaguerreL(q::Complex, p::Complex, x::Complex)
    LaguerreL(CF(real(q), imag(q)), CF(real(p), imag(p)), CF(real(x), imag(x)))
end
function LaguerreL(q::arb, p::arb, x::arb)
    LaguerreL(CF(q), CF(p), CF(x))
end
function LaguerreL(q::Real, p::Real, x::Real)
    LaguerreL(CF(q), CF(p), CF(x))
end
################# Guseinoc Laguerre polynomials L_{p}^{q}(x)
function LaguerreG(q::acb, p::acb, x::acb)
    Power(CF(-1), p) * Gamma(q + CF(1)) * LaguerreL(q-p, p, x)
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                          SPHERICAL HARMONICS                                            #
###########################################################################################################
###########################################################################################################
###########################################################################################################
########################################### SPHERICAL HARMONICS ###########################################(43)
################# Condon-Shortley phases
################# Complex Spherical Harmonics
function SphericalHarmonicsY(n::UInt, m::UInt, theta::acb, phi::acb)
    z = parent(theta)()
    ccall((:acb_hypgeom_spherical_y, Nemo.libarb), Nothing, (Ref{acb}, UInt, UInt, Ref{acb}, Ref{acb}, Int),
    z, n, m, theta, phi, parent(theta).prec)
    return z
end
function SphericalHarmonicsY(n::Int, m::Int, theta::acb, phi::acb)
    if m >= 0
        SphericalHarmonicsY(UInt(n), UInt(m), theta, phi)
    else
        m = abs(m)
        SphericalHarmonicsY(UInt(n), -UInt(m), theta, phi)
    end
end
function SphericalHarmonicsY(n::Int, m::Int, theta::arb, phi::arb)
    SphericalHarmonicsY(n, m, CF(theta), CF(phi))
end
function SphericalHarmonicsY(n::Int, m::Int, theta::Real, phi::Real)
    SphericalHarmonicsY(n, m, CF(theta), CF(phi))
end

function SphericalHarmonicsY(n::UInt, m::UInt, x::arb, y::arb, z::arb)
    r = Sqrt(x^2 + y^2 + z^2)
    SphericalHarmonicsY(n, m, r * Sin(theta) * Cos(phi), r * Sin(theta) * Sin(phi), r * Cos(theta))    
end
################# Complex Spherical Harmonics Guseinov formulae
function SphericalHarmonicsYG(n::Int, m::Int, x::arb, y::arb, z::arb)

    lowl = Zero(0)
    upl = Int(floor((n-abs(m))/2))
    res2 = zeros(RF, upl - lowl + 1)
    r = Sqrt(x^2 + y^2 + z^2)

    coef = Sqrt((2*n + 1) // (4*Pi(RF)))
    res1 = Power((x // r) + CF(0,1)* Epsilon(m,0) * ( y // r), Abs(m))

    for k in lowl : upl
        res2[k - lowl + 1] = SphCB(k, n, abs(m))*((z // r)^(n-abs(m)-2*k))
    end
    coef * res1 * sum(res2)

end

function SphericalHarmonicsYG(n::Int, m::Int, theta::arb, phi::arb)
    lowl = Zero(0)
    upl = Int(floor((n-abs(m))/2))
    res2 = zeros(RF, upl - lowl + 1)

    coef = Sqrt((2*n + 1) // (4*Pi(RF)))
    res1 = Power((Sin(theta) * Cos(phi)) + CF(0,1)* Epsilon(m,0) * (Sin(theta) * Sin(phi)), Abs(m))

    for k in lowl : upl
        res2[k - lowl + 1] = SphCB(k, n, abs(m))*((Cos(theta))^(n-abs(m)-2*k))
    end
    coef * res1 * sum(res2)
end
###########################################################################################################
################# SphericalHarmonicsS already in Condon-Shortley phase
function SphericalHarmonicsS(n::Int, m::Int, theta::acb, phi::acb)
    res1 = RSqrt(Pi(CF) * (One(m) + KroneckerDelta(m, Zero(m))))
    res2 = NLegendreP(CF(n), Abs(CF(m)), Cos(theta), UInt(0))
    res3 = ifelse(m >= 0 , Cos(Abs(m)*phi), Sin(Abs(m)*phi))
    res1 * res2 * res3
end
function SphericalHarmonicsS(n::Int, m::Int, theta::arb, phi::arb)
    SphericalHarmonicsS(n, m, CF(theta), CF(phi))
end
function SphericalHarmonicsS(n::Int, m::Int, theta::Real, phi::Real)
    SphericalHarmonicsS(n, m, CF(theta), CF(phi))
end

function SphericalHarmonicsS(n::Int, m::Int, x::arb, y::arb, z::arb)
    r = Sqrt(x^2 + y^2 + z^2)
    theta = ArcTan(y // x)
    phi = ArcCos(z // r)
    SphericalHarmonicsS(n, m, theta, phi)  
end
################# Real Spherical Harmonics Guseinov formulae
################# ((-1)^m) for Condon-Shortley phase
function SphericalHarmonicsSG(n::Int, m::Int, theta::arb, phi::arb)

    lowl = Zero(0)
    upl = Int(floor((n-abs(m))/2))
    res2 = zeros(RF, upl - lowl + 1)

    if m > 0
        coef = (1 // 2) * Sqrt((2*n + 1) // (2 * Pi(RF)))
    elseif m < 0
        coef = -CF(0,1) * (1 // 2) * Sqrt((2*n + 1) // (2 * Pi(RF)))
    else
        coef = (1 // 2) * Sqrt((2*n + 1) // (1 * Pi(RF)))
    end

    xdr = Sin(theta) * Cos(phi)
    ydr = Sin(theta) * Sin(phi)
    zdr = Cos(theta)
    resxyp = Power(xdr + CF(0,1) * (ydr), Abs(m))
    resxyn = Power(xdr - CF(0,1) * (ydr), Abs(m))

    if m > 0
        res1 = resxyp + resxyn
    elseif m < 0
        res1 = resxyp - resxyn
    else
        res1 = RF(1)
    end

    for k in lowl : upl
        res2[k - lowl + 1] = SphCB(k, n, abs(m))*((zdr)^(n-abs(m)-2*k))
    end

    coef * res1 * sum(res2)
end

function SphericalHarmonicsSG(n::Int, m::Int, x::arb, y::arb, z::arb)

    lowl = Zero(0)
    upl = Int(floor((n-abs(m))/2))
    res2 = zeros(RF, upl - lowl + 1)

    if m > 0
        coef = (1 // 2) * Sqrt((2*n + 1) // (2 * Pi(RF)))
    elseif m < 0
        coef = -CF(0,1) * (1 // 2) * Sqrt((2*n + 1) // (2 * Pi(RF)))
    else
        coef = (1 // 2) * Sqrt((2*n + 1) // (1 * Pi(RF)))
    end

    r = Sqrt(x^2 + y^2 + z^2)
    xdr = x // r
    ydr = y // r
    zdr = z // r
    resxyp = Power(xdr + CF(0,1) * (ydr), Abs(m))
    resxyn = Power(xdr - CF(0,1) * (ydr), Abs(m))

    if m > 0
        res1 = resxyp + resxyn
    elseif m < 0
        res1 = resxyp - resxyn
    else
        res1 = RF(1)
    end

    for k in lowl : upl
        res2[k - lowl + 1] = SphCB(k, n, abs(m))*((zdr)^(n-abs(m)-2*k))
    end
    coef * res1 * sum(res2)

end

function SphericalHarmonicsSG(n::Int, m::Int, x::Real, y::Real, z::Real)
    SphericalHarmonicsSG(n, m, RF(x), RF(y), RF(z))
end
###########################################################################################################
############################################ SPHERICAL SPINORS ############################################(43)
function SphericalSpinors(β::Int, κ::Int, m::arb, theta::arb, phi::arb)
    lᵦ = ifelse(β*κ > 0, β*κ - 1, -β*κ)
    lᵦ = Int(convert(Float64, lᵦ))
    ml0 = Int(convert(Float64,mλ(m,0)))
    ml1 = Int(convert(Float64,mλ(m,1)))
    ucg = Sign(β * κ) * Sqrt((-β*κ + 1//2 - m) // (-2*β*κ + 1))
    lcg = Sqrt((-β*κ + 1//2 + m) // (-2*β*κ + 1))
    uc = ucg * Eta(mλ(m, 0)) * SphericalHarmonicsY(lᵦ, ml0, theta, phi)
    lc = lcg * Eta(mλ(m, 1)) * SphericalHarmonicsY(lᵦ, ml1, theta, phi)
    [uc, lc]
end

function SphericalSpinors(β::Int, κ::Int, m::arb, x::arb, y::arb, z::arb)
    r = Sqrt(x^2 + y^2 + z^2)
    theta = ArcTan(y // x)
    phi = ArcCos(z // r)

    lᵦ = ifelse(β*κ > 0, β*κ - 1, -β*κ)
    lᵦ = Int(convert(Float64, lᵦ))
    ml0 = Int(convert(Float64,mλ(m,0)))
    ml1 = Int(convert(Float64,mλ(m,1)))
    ucg = Sign(β * κ) * Sqrt((-β*κ + 1//2 - m) // (-2*β*κ + 1))
    lcg = Sqrt((-β*κ + 1//2 + m) // (-2*β*κ + 1))
    uc = ucg * Eta(mλ(m, 0)) * SphericalHarmonicsY(lᵦ, ml0, theta, phi)
    lc = lcg * Eta(mλ(m, 1)) * SphericalHarmonicsY(lᵦ, ml1, theta, phi)
    [uc, lc]
end

function RadoslawSSpinors(κ::Int, m::arb, theta::arb, phi::arb)

    l = ifelse(κ > 0, κ, -κ - 1)
    ml0 = Int(convert(Float64,mλ(m,0)))
    ml1 = Int(convert(Float64,mλ(m,1)))
    ucg = Sign(-κ) * Sqrt((κ + 1//2 -m) // (2 * κ + 1))
    lcg = Sqrt((κ + 1//2 + m) // (2 * κ + 1))
    uc = ucg * Eta(mλ(m, 0)) * SphericalHarmonicsY(l, ml0, theta, phi)
    lc = lcg * Eta(mλ(m, 1)) * SphericalHarmonicsY(l, ml1, theta, phi)
    [uc, lc]
end

function SphericalSpinors(l::Int, j::arb, m::arb, theta::arb, phi::arb)
    κ = 2*(l-j)*(j + 1/2)
    κ = Int(convert(Float64, κ))
    uc = SphericalSpinors(-1, κ, m, theta, phi)
    lc = SphericalSpinors(1, κ, m, theta, phi)
    reshape([uc lc], 4,1)
end

function SphericalSpinors(l::Int, j::arb, m::arb, x::arb, y::arb, z::arb)
    κ = 2*(l-j)*(j + 1/2)
    κ = Int(convert(Float64, κ))
    uc = SphericalSpinors(-1, κ, m, x, y, z)
    lc = SphericalSpinors(1, κ, m, x, y, z)
    reshape([uc lc], 4,1)
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                     ROTATED ANGULAR FUNCTIONS                                           #
###########################################################################################################
###########################################################################################################
###########################################################################################################(43)
################# x = 0 for the complex spherical harmonics, x = 1 for the real ones
function RotaD(λ::Int, l1::Int, m1::Int, l2::Int, m2::Int, theta::arb, phi::arb)
    coef = (RF(2) // (RF(1) + KroneckerDelta(λ,0)))
    resr = RF(0)
    for L in abs(l1-l2) : 2 : l1 + l2 
        res1 = ClebschGordanG(l1, -m1, l2, m2, L, (-m1 + m2))
        res2 = ClebschGordanG(l1, -λ, l2, λ, L, 0)
        res3 = Sqrt((RF(4) * Pi(RF)) // (2 * L + 1))
        res4 = Eta(-m1 + m2) * SphericalHarmonicsY(L, -m1 + m2, theta, phi)
        resr += res1 * res2 * res3 * res4
    end
    res = coef * resr
    return res
end

function Rotad(λ::Int, l1::Int, m1::Int, l2::Int, m2::Int, theta::arb, phi::arb)
    if (theta == 0) && (phi == 0)
        res = KroneckerDelta(m1,m2) * KroneckerDelta(λ, abs(m1))
    elseif (theta == Pi(RF)) && (phi == 0)
        res = KroneckerDelta(m1,m2) * KroneckerDelta(λ, abs(m1)) * ((-1)^(l1+l2))
    else
        coef1 = (
            Power(CF(0,1), m1 >= 0 ? 1 : 0) *
            Power(CF(0,-1),m2 >= 0 ? 1 : 0)
        )
        coef2 = Sqrt(
            (RF(1) + KroneckerDelta(m1,0)) *
            (RF(1) + KroneckerDelta(m2,0))
        )
        coef = coef1//(2*coef2)
        RotaD1 = RotaD(λ, l1, abs(m1), l2, abs(m2), theta, phi)
        RotaD2 = Epsilon(m1, 0) * RotaD(λ, l1, -abs(m1), l2, abs(m2), theta, phi)
        RotaD3 = Epsilon(0, m2) * RotaD(λ, l1, abs(m1), l2, -abs(m2), theta, phi)
        RotaD4 = Epsilon(m1, m2) * RotaD(λ, l1, -abs(m1), l2, -abs(m2), theta, phi)
        res = coef * (RotaD1 + RotaD2 + RotaD3 + RotaD4)
    end
    return res
end

function RotaT(λ::Int, l1::Int, m1::Int, l2::Int, m2::Int, theta::arb, phi::arb, x::Int)
    if x == 0
        RotaD(λ::Int, l1::Int, m1::Int, l2::Int, m2::Int, theta::arb, phi::arb)
    else
        Rotad(λ::Int, l1::Int, m1::Int, l2::Int, m2::Int, theta::arb, phi::arb)
    end
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                          SLATER-TYPE FUNCTIONS                                          #
###########################################################################################################
###########################################################################################################
###########################################################################################################
function SlaterORadial(n::arb, ζ::arb, r::arb)
    (r^(n - 1)) * Exp(-ζ * r)
end
function SlaterORadial(n::Real, ζ::Real, r::Real)
    SlaterORadial(RF(n), RF(ζ), RF(r))
end

function SlaterOrbitalY(n, l, m, ζ, r, theta, phi)
    SlaterON(n, ζ) * SlaterORadial(n, ζ, r) * SphericalHarmonicsY(l, m, theta, phi)
end
function SlaterOrbitalY(n, l, m, ζ, r, x, y, z)
    r = Sqrt(x^2 + y^2 + z^2)
    theta = ArcTan(y // x)
    phi = ArcCos(z // r)

    SlaterON(n, ζ) * SlaterORadial(n, ζ, r) * SphericalHarmonicsY(l, m, theta, phi)
end

function SlaterOrbitalS(n, l, m, ζ, r, theta, phi)
    SlaterON(n, ζ) * SlaterORadial(n, ζ, r) * SphericalHarmonicsS(l, m, theta, phi)
end
function SlaterOrbitalS(n, l, m, ζ, r, x, y, z)
    r = Sqrt(x^2 + y^2 + z^2)
    theta = ArcTan(y // x)
    phi = ArcCos(z // r)

    SlaterON(n, ζ) * SlaterORadial(n, ζ, r) * SphericalHarmonicsS(l, m, theta, phi)
end
###########################################################################################################
###################################### CHARGE DENSITY DISTRIBUTION ########################################
function ChargeDOneC(
    n1::arb, l1::Int, m1::Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    N::arb, L::Int, M::Int
    )
    
    z = ζ1 + ζ2
    τ = ((ζ1 - ζ2) // (ζ1 + ζ2))

    res1 = (Power(z, 3//2) // Power(2, N))
    res2 =  ((2*L +1) // 2)
    res3 = ((Gamma(2*N + 1)) // (Gamma(2*n1 + 1) * Gamma(2*n2 + 1)))
    res4 = (1 + τ)
    res5 = (1 - τ)
    res6 = GGauntG(l1, m1, l2, m2, L, abs(M))
    res7 = RSphCA(M, m1, m2)

    res = res1 * Sqrt(res2 * res3) * res4 * res5 * res6 * res7
    return res
end
###########################################################################################################
######################################### ONE-CENTER POTENTIAL ############################################
function OneCenterP(L::Int, n1::arb, n2::arb, ζ1::arb, ζ2::arb, R::arb)
    ρ = (R//2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)

    coeff1 = SlaterON(n1, n2, RF(1), τ) * (2 * (ζ1 + ζ2))
    coeff2 = Gamma(n1 + n2 + L + 1)
    coeff3 = Power(2 * ρ, L + 1)
    coeff = coeff1 * coeff2 * (RF(1) // coeff3)

    res1 = PGamma(n1 + n2 + L + 1, 2 * ρ)
    res2 = Power(2 * ρ, 2 * L + 1)
    res3 = Pochhammer(n1 + n2 - L, 2 * L + 1)
    res4 = QGamma(n1 + n2 - L, 2 * ρ)
    res5 = res1 + (res2 // res3) * res4

    res = coeff * res5

    return res
end
###########################################################################################################
########################################### SLATER TYPE SPINORS ###########################################
################# Two-component form radial
function SlaterSRadial(β::Int, n::arb, κ::Int, ζ::arb, r::arb)
    ((SlaterSA(β, n, κ) * (r^n)) + (SlaterSB(β, n, κ) * (r^(n-1)))) * Exp(-ζ * r)
end

function SlaterSRadial(β::Int, n::arb, l::Int, j::arb, ζ::arb, r::arb)
    κ = 2*(l-j)*(j + 1/2)
    κ = Int(convert(Float64, κ))
    ((SlaterSA(β, n, κ) * (r^n)) + (SlaterSB(β, n, κ) * (r^(n-1)))) * Exp(-ζ * r)
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                          HYPERGEOMETRIC FUNCTIONS                                       #
###########################################################################################################
###########################################################################################################
############################################ HYPERGEOMETRIC0F0 ############################################(44)
Hypergeometric0F0(x::acb) = Exp(x)
Hypergeometric0F0(x::Complex) = Exp(x)
Hypergeometric0F0(x::arb) = Exp(x)
Hypergeometric0F0(x::Real) = Exp(x)
###########################################################################################################
############################################ HYPERGEOMETRIC1F0 ############################################(44)
Hypergeometric1F0(a::acb, x::acb) = Surd(1-x, a)
Hypergeometric1F0(a::Complex, x::Complex) = Surd(1-x, a)
Hypergeometric1F0(a::arb, x::arb) = Surd(1-x, a)
Hypergeometric1F0(a::Real, x::Real) = Surd(1-x, a)
###########################################################################################################
############################################ HYPERGEOMETRIC0F1 ############################################(44)
################# Confluent hypergeometric limit functions
################# Confluent hypergeometric function 0f1. For r = (0,1) = (irregularized, regularized)
################# r = 1 -> 0F1 = 0F1/Gamma(a)
function Hypergeometric0F1(a::acb, x::acb, r::Int)
    z = parent(x)()
    ccall((:acb_hypgeom_0f1, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
    z, a, x, r, parent(x).prec)
    return z
end
Hypergeometric0F1(a::Complex, x::Complex, r::Int) =
Hypergeometric0F1(CF(real(a), imag(a)), CF(real(x), imag(x)), r)

Hypergeometric0F1(a::arb, x::arb, r::Int) = Hypergeometric0F1(CF(a), CF(x), r)
Hypergeometric0F1(a::Real, x::Real, r::Int) = Hypergeometric0F1(CF(a), CF(x), r)
###########################################################################################################
############################################ HYPERGEOMETRIC1F1 ############################################(44)
################# Confluent hypergeometric functions of the first kind also writtten as M(a;b;x)
################# Kummer confluent hypergeometric functions
################# (0,1) = (irregularized, regularized)
function Hypergeometric1F1(a::acb, b::acb, x::acb, r::Int)
    z = parent(x)()
    ccall((:acb_hypgeom_m, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
    z, a, b, x, r, parent(x).prec)
    return z
end
Hypergeometric1F1(a::Complex, b::Complex, x::Complex, r::Int) =
Hypergeometric1F1(CF(real(a), imag(a)), CF(real(b), imag(b)), CF(real(x), imag(x)), r)

Hypergeometric1F1(a::arb, b::arb, x::arb, r::Int) = Hypergeometric1F1(CF(a), CF(b), CF(x), r)
Hypergeometric1F1(a::Real, b::Real, x::Real, r::Int) = Hypergeometric1F1(CF(a), CF(b), CF(x), r)
###########################################################################################################
############################################# HYPERGEOMETRICU #############################################(45)
function HypergeometricU(a::acb, b::acb, x::acb)
    z = parent(x)()
    ccall((:acb_hypgeom_u, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Int),
    z, a, b, x, parent(x).prec)
    return z
end
HypergeometricU(a::Complex, b::Complex, x::Complex) =
HypergeometricU(CF(real(a), imag(a)), CF(real(b), imag(b)), CF(real(x), imag(x)))

HypergeometricU(a::arb, b::arb, x::arb) = HypergeometricU(CF(a), CF(b), CF(x))
HypergeometricU(a::Real, b::Real, x::Real) = HypergeometricU(CF(a), CF(b), CF(x))
###########################################################################################################
################################################ WHITTAKER ################################################(48)
################# WhittakerM(k,m,z) = Exp(-z/2)*Power(z, m+1/2)*1F1(1/2+m-k,1+2m;z)
function WhittakerM(k, m, x) 
    Exp(-x//2) * Power(x, m + (1//2)) * Hypergeometric1F1(1//2 + m - k, One(m) + 2*m, x, sprec)
end
################# WhittakerW(k,m,z) = Exp(-z/2)*Power(z, m+1/2)*U(1/2+m-k,1+2m;z)
function WhittakerW(k, m, x)
    Exp(-x//2) * Power(x, m + (1//2)) * HypergeometricU(1//2 + m - k, One(m) + 2*m, x)
end
###########################################################################################################
############################################ HYPERGEOMETRIC2F1 ############################################(44)
function Hypergeometric2F1(a::acb, b::acb, c::acb, x::acb, r::Int)
    z = parent(x)()
    ccall((:acb_hypgeom_2f1, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
    z, a, b, c, x, r, parent(x).prec)
    return z
end
Hypergeometric2F1(a::Complex, b::Complex, c::Complex, x::Complex, r::Int) = Hypergeometric2F1(
    CF(real(a), imag(a)), 
    CF(real(b), imag(b)), 
    CF(real(c), imag(c)), 
    CF(real(x), imag(x)), 
    r
    )
Hypergeometric2F1(a::arb, b::arb, c::arb, x::arb, r::Int) = Hypergeometric2F1(
    CF(a), CF(b), CF(c), CF(x), r
    )
Hypergeometric2F1(a::Real, b::Real, c::Real, x::Real, r::Int) = Hypergeometric2F1(
    CF(a), CF(b), CF(c), CF(x), r
    )
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                     MOLECULAR AUXILIARY FUNCTIONS                                       #
###########################################################################################################
###########################################################################################################
############################################## OVERLAP-LİKE G #############################################(45)
################# Overlap-like G used for calculationg the two-center one-electron integrals
###########################################################################################################
###########################################################################################################
################# Case 1 The parameter p_{3} = 0
function AuxiliaryGk(n1::arb, q::Int, n2::arb, n3::arb, p1::arb, p2::arb)
    two = parent(n1)(2)
    res1 = ((p1^n1) // (Gamma(n1+1))) * (two ^ (n2+n3+q+One(n2)))
    res2 = Beta(n2+1, n3+1, parent(n2)(1/2)) * ExpIntegralE(-(n2+n3+q+1), p2)
    res = res1 * res2
end

function AuxiliaryGm(n1::arb, n2::arb, p::arb)
    two = parent(n1)(2)
    (two ^ n1) * HypergeometricU(n2+One(n2), n1+n2+two, p) * Gamma(n2+One(n2)) * Exp(-p)
end

function AuxiliaryGl(n1::arb, q::Int, n2::arb, n3::arb, p1::arb, p2::arb, lim::Int)
    two = parent(n1)(2)
    lowl = Zero(lim)
    upl = lim

    a = (n3) // (n2 + q + One(n2))
    b = (n2-n3+q+One(n2)-p2) // (n2+q+One(n2))

    res = zeros(RF, upl-lowl+1)
    p = zeros(RF, upl-lowl+1)
    f = zeros(RF, upl-lowl+1)
    m1 = zeros(RF, upl-lowl+1)
    m2 = zeros(RF, upl-lowl+1)
    m = zeros(RF, upl-lowl+1)

    p[1] = Pochhammer(-n2,0)
    f[1] = RF(1)
    m1[1] = real(AuxiliaryGm(n2 + q + two, n3 - One(n3), p2))
    m2[1] = real(AuxiliaryGm(n2 + q + One(n2), n3, p2))
    m[1] = (1//4) * a * m1[1] + (1//2) * b * m2[1]

    for s in lowl : upl-1
        s1 = s-lowl+1
        p[s1+1] = (-n2 + s1 -One(n2)) * p[s1]
        f[s1+1] = (s+1) * f[s1]
        a = (n3+s1) // (n2 + q + One(n2) - s1)
        b = (n2-n3-2*s1+q+One(n2)-p2) // (n2+q+One(n2)-s1)
        m[s1+1] = (1//4) * a * m2[s1] + (1//2) * b * m[s1]
        m2[s1+1] = m[s1]
    end
    for s in lowl : upl
        s1 = s-lowl+1
        res[s1] = (p[s1] // ((n3 + s + One(n3)) * f[s1])) * m[s1]
    end
    ((p1^n1) // Gamma(n1 + One(n1))) * sum(res)
end

function AuxiliaryGh(n1::arb, q::Int, n2::arb, n3::arb, p1::arb, p2::arb, lim::Int)
    two = parent(n1)(2)
    res1 = ((p1^n1) // Gamma(n1 + One(n1))) * (two^(n2+n3+q+One(n1))) * Beta(n2+One(n2), n3+One(n3))
    res2 = ExpIntegralE(-(n2 + n3 + q + One(n1)), p2)
    res3 = AuxiliaryGl(n1, q, n2, n3, p1, p2, lim)
    res1 * res2 - res3
end

function AuxiliaryG(n1::arb, n2::arb, n3::arb, p1::arb, p2::arb, lim::Int)
    two = parent(n1)(2)
    res1 = ((p1^n1) // Gamma(n1 + One(n1))) * (two^(n2+n3+One(n1))) * Beta(n2+One(n2), n3+One(n3))
    res2 = ExpIntegralE(-(n2 + n3 + One(n1)), p2)
    res3 = AuxiliaryGl(n1, 0, n2, n3, p1, p2, lim)
    res4 = AuxiliaryGl(n1, 0, n3, n2, p1, p2, lim)
    res1 * res2 - res3 - res4
end

function AuxiliaryGrlm(n1::arb, n2::arb, n3::arb, p1::arb, p2::arb, p3::arb)
    two = parent(n1)(2)
    res1 = ((p1^n1) // Gamma(n1 + One(n1))) * (two^(n3 + One(n1))) * Exp(-p2) * Exp(p3)
    res2 = AuxiliaryGm(n2, n3, two*p3)
    res1*res2
end

function AuxiliaryGrlp(n1::arb, n2::arb, n3::arb, p1::arb, p2::arb, p3::arb)
    two = parent(n1)(2)
    res1 = ((p1^n1) // Gamma(n1 + One(n1))) * (two^(n3+One(n1))) * Exp(p2) * Exp(p3)
    res2 = AuxiliaryGm(n2, n3, two*p3)
    res1*res2
end

function AuxiliaryGrh(n1::arb, n2::arb, n3::arb, p1::arb, p2::arb, p3::arb, prec::Int)
    two = parent(n1)(2)
    res1 = ((p1^n1) // Gamma(n1 + One(n1))) * (two^(n2 + n3 + One(n2 + n3)))
    res2 = Gamma(n2 + One(n2)) * Gamma(n3 + One(n3))
    res3 = Hypergeometric1F1(n2 + One(n2), n2 + n3 + two, -two * p3, prec)
    res = Exp(p3 - p2) * res1 * res2 * res3
end
###########################################################################################################
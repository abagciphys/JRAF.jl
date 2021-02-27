###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                           ANGULAR COEFFICIENTS                                          #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
export Epsilon, Eta, mλ, RSphCA, CSphCA, SphCB, SphPCD, SphPCG0, SphPCG, ClebschGordanG, ClebschGordan, 
GauntGF, GauntG, GGauntG, GauntWF, GauntW
###########################################################################################################
############################################ SIGN COEFFICIENTS ############################################
Epsilon(x, y) = Times(Nonnegative(x), Nonnegative(y))
Eta(x) = Power((-1), ((abs(x)-x)//2)) # Spherical spinors sign coefficient to convert to Bethe phase
Eta(x, y, z) = Sign(x, y, z)
mλ(m, λ) = m + λ - 1//2 # Nonrelativistic angular momentum quamtum number m from mj
###########################################################################################################
############################### SPHERICAL HARMONICS PRODUCT COEFFICIENTS A ################################(25)
################# A_{σ}(m, mü)
function RSphCA(σ::Int, m1::Int, m2::Int)
    res1 = Sqrt(2-Abs(Eta(m1, m2, m1-m2)))*KroneckerDelta(σ, Epsilon(m1, m2)*Abs(m1-m2))
    res2 = Eta(m1, m2, m1+m2)*KroneckerDelta(σ, Epsilon(m1, m2)*Abs(m1+m2))
    RSqrt(2) * (res1 + res2)
end
RSphCA(σ::arb, m1::arb, m2::arb) = RSphCA(
    Int(convert(Float64, σ)),
    Int(convert(Float64, m1)),
    Int(convert(Float64, m2))
)
RSphCA(σ::Real, m1::Real, m2::Real) = RSphCA(Int(σ), Int(m1), Int(m2))

CSphCA(σ::Int, m1::Int, m2::Int) = σ == m1 - m2 ? One(ArbField(sprec)) : Zero(ArbField(sprec))
CSphCA(σ::arb, m1::arb, m2::arb) = CSphCA(
    Int(convert(Float64, σ)),
    Int(convert(Float64, m1)),
    Int(convert(Float64, m2))
)
CSphCA(σ::Real, m1::Real, m2::Real) = CSphCA(Int(σ), Int(m1), Int(m2))
###########################################################################################################
################################### SPHERICAL HARMONICS COEFFICIENTS B ####################################(25)
################# B(k,l,m)
function SphCB(k::Int, l::Int, λ::Int)
    km2 = 2 * k
    lm2 = 2 * l

    res1 = (((-1)^k)//Power(2,l))*RSqrt(Binomial(l, λ) * Binomial(l + λ, λ))
    res2 = Binomial(k + λ, k) * Binomial(lm2 - km2, l - k) * Binomial(l - k, l - λ - km2)
    res1 * res2
end
SphCB(k::arb, l::arb, m::arb) = SphCB(
    Int(convert(Float64, k)),
    Int(convert(Float64, l)),
    Int(convert(Float64, m))
)
SphCB(k::Real, l::Real, λ::Real) = SphCB(Int(k), Int(l), Int(λ))
###########################################################################################################
############################### SPHERICAL HARMONICS PRODUCT COEFFICIENTS D ################################(25)
################# D(l,lambda,beta)
function SphPCD(l::Int, λ::Int, β::Int)
    lmβd2 = (l - β) / 2 # Exponentiation below line yielding a complex result when mod(l - b, 2) != 0
    #res1 = ((-1)^lmβd2) / (2^l)
    res1 = (EvenQ(lmβd2) ? 1 : -1) / (2^l)
    res2 = Sqrt((2 * l + 1) * Binomial(l + λ, l) / (2 * Binomial(l, λ)))
    res3 = Binomial(l, floor(Int, lmβd2)) * Binomial(l + β, β - λ)
    res1 * res2 * res3
end
SphPCD(l::arb, λ::arb, β::arb) = β < 0 ? 0 : SphPCD(
    UInt(convert(Float64, l)),
    UInt(convert(Float64, λ)),
    UInt(convert(Float64, β))
)
SphPCD(l::Real, λ::Real, β::Real) = β < 0 ? 0 : SphPCD(Int(l), Int(λ), Int(β))
###########################################################################################################
############################### SPHERICAL HARMONICS PRODUCT COEFFICIENTS G ################################(25)
################# G(alpha, beta, l1, lambda1, l2, lambda2, Lambda)
function SphPCG0(α::Int, β::Int, l1::Int, λ1::Int, l2::Int, λ2::Int, Λ::Int)
    res = zeros(RF, Λ + 1)
    for s in 0 : Λ
        res[s + 1] = ((-1)^(s)) * Binomial(Λ, s) * SphPCD(l1, λ1, α + 2 * Λ - 2 * s)
    end
    sum(res) * SphPCD(l2, λ2, β)
end
#################
function SphPCG(q::Int, α::Int, β::Int, l1::Int, λ1::Int, l2::Int, λ2::Int, Λ::Int)
    SphPCG0(α, β, l1, λ1, l2, λ2, Λ) * GBinomial(q, α + 2 * Λ - λ1, β - λ2)
end
###########################################################################################################
####################################### CLEBSCH-GORDAN COEFFICIENTS #######################################(32)
################# I. I. Guseinov formulation
function ClebschGordanG(l1::Int, m1::Int, l2::Int, m2::Int, L::Int, M::Int)
    if L < M || m1 + m2 != M
        return Zero(ArbField(sprec))
    else
        phs = (-1)^((m1 + abs(m1) + m2 + abs(m2) + M + abs(M)) / 2)
        Lt2 = 2 * L + 1
        l1t2 = 2 * l1 + 1
        l2t2 = 2 * l2 + 1

        lL1 = l1 + l2 + L + 1
        lL2 = l1 + l2 - L
        lL3 = l1 - l2 + L
        lL4 = l2 - l1 + L
        res1 = (Lt2^2)*Binomial(lL1, lL2)*Binomial(2*L, L+M)
        res2 = l1t2*l2t2*Binomial(lL1, lL3)*Binomial(lL1, lL4)*Binomial(2*l1, l1+m1)*Binomial(2*l2, l2+m2)

        lowl = max(0, l1-m1-L+M, l2+m2-L-M)
        upl = min(l1-m1, l2+m2, l1+l2-L)
        
        if lowl > upl
            lowl = upl = 0
        end

        res3 = zeros(RF, upl-lowl+1)
        a1 = zeros(RF, upl-lowl+1)
        a2 = zeros(RF, upl-lowl+1)
        a3 = zeros(RF, upl-lowl+1)
        a1[1] = Binomial(lL2,lowl)
        a2[1] = Binomial(L-M,l1-m1-lowl)
        a3[1] = Binomial(L+M,l2+m2-lowl)

        for s in lowl : upl-1
            s1 = s-lowl+1
            a1[s1+1] = ((lL2-s) // (s+1))*a1[s1]
            a2[s1+1] = ((l1-m1-s) // ((L-M)-(l1-m1-s)+1))*a2[s1]
            a3[s1+1] = ((l2+m2-s) // ((L+M)-(l2+m2-s)+1))*a3[s1]
        end

        for s in lowl : upl
            s1 = s-lowl+1
            res3[s1] = ((-1)^s)*a1[s1]*a2[s1]*a3[s1]
        end
        res = phs * KroneckerDelta(M,m1 + m2) * rsqrt(res2) * sqrt(res1) * sum(res3)
        return res
    end
end
################# Using Wigner3j symbol via WignerSymbol package
function ClebschGordan(l1, m1, l2, m2, L, M)
    return clebschgordan(BigFloat, l1, m1, l2, m2, L, M)
end
###########################################################################################################
########################################### GAUNT COEFFICIENTS ############################################(37)
################# If "," is used the given value is defauls. If ";" is used the value assigned.
function GauntGF(l1::Int, m1::Int, l2::Int, m2::Int, L::Int; M::Int = m1-m2)
    lowl = max(0, L-m1-l2)
    upl = min(l1-abs(m1), L-M, L-m1+l2)

    if EvenQ(l1 + l2 + L) == true && l1 - l2 <= L <= l1 + l2 && upl >= lowl
        γ = (l1 + l2 + L) // 2

        phs = (-1)^(γ - (l2-m2) + (abs(m1)+abs(m2)+abs(M)) // 2)
        res1 = (Binomial(2*γ-l1-l2, γ-l1)*Binomial(γ, L)) // ((2*γ+1)Binomial(2*γ, 2*L))
        res2 = (2*l1+1)*(2*l2+1)*Binomial(l1+l2+M, l1+m1)*Binomial(2*L+l1+l2+M, l1+l2+M)
        res3 = Binomial(l1+l2-M, l1-m1)*Binomial(2*L+l1+l2+M, l1+l2-M)*Binomial(2*L, L-M)*Binomial(2*L+2*M, L+M)
        res4 = Binomial(l1+l2+M, l1+l2-L) // Binomial(l1+l2+M, l1+m1)

        bin1 = zeros(RF, upl-lowl+1)
        bin2 = zeros(RF, upl-lowl+1)
        bin3 = zeros(RF, upl-lowl+1)
        bin4 = zeros(RF, upl-lowl+1)

        bin1[1] = Binomial(l1+m1+lowl,lowl)
        bin2[1] = Binomial(l1+l2-L, l1-m1-lowl)
        bin3[1] = Binomial(l2-m2+L-M-lowl, L-M-lowl)
        bin4[1] = ((-1)^lowl) * bin1[1] * bin2[1] * bin3[1]

        for s in lowl : upl-1
            s1 = s-lowl+1
            bin1[s1+1] = ((l1+m1+s+1) // (s+1))*bin1[s1]
            bin2[s1+1] = ((l1-m1-s) // (l2-L+m1+s+1))*bin2[s1]
            bin3[s1+1] = ((L-M-s) // (l2-m2+L-M-s))*bin3[s1]
            bin4[s1+1] = ((-1)^(s+1)) * bin1[s1+1] * bin2[s1+1] * bin3[s1+1]
        end
        res5 = res1 * (Sqrt(res2 // res3)) * res4
        res = phs * res5 * sum(bin4)
        return res
    else
        return Zero(ArbField(sprec))
    end
end

function GauntG(l1::Int, m1::Int, l2::Int, m2::Int, L::Int)
    M = m1-m2
    if l1 >= l2
        if m1 >= 0 && m2 >= 0
            GauntGF(l1, m1, l2, m2, L)
        elseif m1 < 0 && m2 >= 0
            GauntGF(l2, m2, l1, m1, L)
        elseif m1 >= 0 && m2 < 0
            Sqrt((2*l1+1) // (2*L+1))*GauntGF(L, M, l2, -m2, l1)
        elseif m1 < 0 && m2 < 0
            GauntGF(l1, -m1, l2, -m2, L)
        else
            Zero(ArbField(sprec))
        end
    else l1 < l2
        if m1 >= 0 && m2 >= 0 || m1 < 0 && m2 >= 0
            GauntGF(l2, m2, l1, m1, L)
        elseif  m1 >= 0 && m2 < 0
            Sqrt((2*l1+1) // (2*L+1))*GauntGF(L, M, l2, -m2, l1)
        elseif m1 < 0 && m2 < 0
            GauntGF(l1, -m1, l2, -m2, L)
        else
            Zero(ArbField(sprec))
        end
    end
end

function GGauntG(l1::Int, m1::Int, l2::Int, m2::Int, L::Int, M::Int)
    if abs(M) == abs(m1 - m2)
        GauntG(l1, m1, l2, m2, L)
    elseif abs(M) == abs(m1 + m2)
        GauntG(l1, m1, l2, -m2, L)
    else
        return Zero(ArbField(sprec))
    end
end
################# Using Wigner3j symbol via WignerSymbol package
function GauntWF(l1::Int, m1::Int, l2::Int, m2::Int, L::Int)
    #phs = (-1)^(((l1 + l2 + L) // 2) - (l2- m1))
    res1 = Sqrt((2*l1 + 1) * (2*l2 + 1)) // (2*L + 1)
    res2 = ClebschGordan(l1, m1, l2, -m2, L, m1 - m2)
    res3 = ClebschGordan(l1, 0, l2, 0, L, 0)
    return res1 * res2 * res3
end

function GauntW(l1::Int, m1::Int, l2::Int, m2::Int, L::Int)
    M = m1 - m2
    if l1 >= l2 && m1 >= 0 && m2 >= 0
        GauntWF(l1, m1, l2, m2, L)
    elseif l1 >= l2 && m1 < 0 && m2 >= 0 || l1 < l2 && m1 >= 0 && m2 >= 0 || l1 < l2 && m1 < 0 && m2 >= 0
        GauntWF(l2, m2, l1, m1, L)
    elseif l1 >= l2 && m1 >= 0 && m2 < 0 || l1 < l2 && m1 >= 0 && m2 < 0
        sqrt((2*l1+1) // (2*L+1))*GauntWF(L, M, l2, -m2, l1)
    elseif l1 >= l2 && m1 < 0 && m2 < 0  || l1 < l2 && m1 < 0 && m2 < 0
        GauntWF(l1, -m1, l2, -m2, L)
    else
        Zero(ArbField(sprec))
    end
end
#For Sage_math results (GGauntG = GauntW) * Sqrt((2 * L + 1) // (4 * Pi(RF)))  
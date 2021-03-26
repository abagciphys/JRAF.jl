###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                         SLATER-TYPE ORBİTALS ONE ELECTRON MOLECULAR INTEGRALS                           #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
export TwoCenterOverlap, CTwoCenterOverlap, VTwoCenterOverlap, STwoCenterOverlap, TwoCenterOverlapRec, 
TwoCenterNucAttractABA, CTwoCenterNucAttractABA, VTwoCenterNucAttractABA, STwoCenterNucAttractABA,
TwoCenterNucAttractAAB, CTwoCenterNucAttractAAB, VTwoCenterNucAttractAAB, VTwoCenterNucAttractAAB,
TwoCenterKinEnergy, CTwoCenterKinEnergy, VTwoCenterKinEnergy, STwoCenterKinEnergy
###########################################################################################################
###########################################################################################################
####################################### TWO CENTER OVERLAP INTEGRALS ######################################(43)
################################################### ALIGNED MOLECULAR COORDINATE SYSTEM
function TwoCenterOverlap(n1::arb, l1::Int, n2::arb, l2::Int, λ::Int, ρ::arb, τ::arb,lim::Int)
    Overlap = RF(0)
    for s1 in (EvenQ(l1-λ) ? -λ : -λ + 1)  : 2 : l1
        for s2 in (EvenQ(l2+λ) ? λ : λ + 1)  : 2 : l2
            for s3 in 0 : s1 + s2
                Overlap += SphPCG(s3, s1, s2, l1, λ, l2, λ, λ)*
                AuxiliaryG(RF(10//10), s3, RF(n1 - s1), RF(n2 - s2), RF(10//10), ρ, ρ * τ, lim)
            end
        end
    end
    return Power(-1,l2 + λ) * SlaterON(n1, n2, ρ, τ) * Overlap
end

function CTwoCenterOverlap(n1::arb, l1::Int, n2::arb, l2::Int, λ::Int, ρ::arb, τ::arb,lim::Float64)
    Overlap = RF(0)
    for s1 in (EvenQ(l1-λ) ? -λ : -λ + 1)  : 2 : l1
        for s2 in (EvenQ(l2+λ) ? λ : λ + 1)  : 2 : l2
            for s3 in 0 : s1 + s2
                Overlap += SphPCG(s3, s1, s2, l1, λ, l2, λ, λ)*
                CuhreAuxiliaryG(10//10, s3, n1 - s1, n2 - s2, 10//10, ρ, ρ * τ, lim)[1]
            end
        end
    end
    return Power(-1,l2 + λ) * SlaterON(n1, n2, ρ, τ) * Overlap
end

function VTwoCenterOverlap(n1::arb, l1::Int, n2::arb, l2::Int, λ::Int, ρ::arb, τ::arb,lim::Float64)
    Overlap = RF(0)
    for s1 in (EvenQ(l1-λ) ? -λ : -λ + 1)  : 2 : l1
        for s2 in (EvenQ(l2+λ) ? λ : λ + 1)  : 2 : l2
            for s3 in 0 : s1 + s2
                Overlap += SphPCG(s3, s1, s2, l1, λ, l2, λ, λ)*
                VegasAuxiliaryG(10//10, s3, n1 - s1, n2 - s2, 10//10, ρ, ρ * τ, lim)[1]
            end
        end
    end
    return Power(-1,l2 + λ) * SlaterON(n1, n2, ρ, τ) * Overlap
end

function STwoCenterOverlap(n1::arb, l1::Int, n2::arb, l2::Int, λ::Int, ρ::arb, τ::arb,lim::Float64)
    Overlap = RF(0)
    for s1 in (EvenQ(l1-λ) ? -λ : -λ + 1)  : 2 : l1
        for s2 in (EvenQ(l2+λ) ? λ : λ + 1)  : 2 : l2
            for s3 in 0 : s1 + s2
                Overlap += SphPCG(s3, s1, s2, l1, λ, l2, λ, λ)*
                SuaveAuxiliaryG(10//10, s3, n1 - s1, n2 - s2, 10//10, ρ, ρ * τ, lim)[1]
            end
        end
    end
    return Power(-1,l2 + λ) * SlaterON(n1, n2, ρ, τ) * Overlap
end

function TwoCenterOverlapRec(n1::arb, l1::Int, n2::arb, l2::Int, λ::Int, ρ::arb, τ::arb,lim::Int)
    OverlapRec = RF(0)
    for s1 in (EvenQ(l1-λ) ? -λ : -λ + 1)  : 2 : l1
        for s2 in (EvenQ(l2+λ) ? λ : λ + 1)  : 2 : l2
            for s3 in 0 : s1 + s2
                OverlapRec += SphPCG(s3, s1, s2, l1, λ, l2, λ, λ)*
                AuxiliaryGr(RF(10//10), s3, RF(n1 - s1), RF(n2 - s2), RF(10//10), ρ, ρ * τ, lim)
            end
        end
    end
    return SlaterON(n1, n2, ρ, τ)*OverlapRec
end
###########################################################################################################
function CTwoCenterOverlap(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64)

    CTwoCenterOverlap(n1, l1, n2, l2, λ, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), lim)
end

function VTwoCenterOverlap(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64)

    VTwoCenterOverlap(n1, l1, n2, l2, λ, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), lim)
end

function STwoCenterOverlap(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64)

    STwoCenterOverlap(n1, l1, n2, l2, λ, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), lim)
end

function TwoCenterOverlapRec(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Int)

    TwoCenterOverlapRec(n1, l1, n2, l2, λ, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), lim)
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################### NON-ALIGNED MOLECULAR COORDINATE SYSTEM
function  TwoCenterOverlap(
n1::arb, l1::Int, m1::Int,
n2::arb, l2::Int, m2::Int,
ρ::arb, τ::arb, theta::arb, phi::arb,
x::Int, lim::Int)

    NOverlap = RF(0)

    for s in 0 : min(l1,l2)
        NOverlap += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        TwoCenterOverlap(n1, l1, n2, l2, s, ρ, τ,lim))
    end
    NOverlap
end

function  CTwoCenterOverlap(
    n1::arb, l1::Int, m1::Int,
    n2::arb, l2::Int, m2::Int,
    ρ::arb, τ::arb, theta::arb, phi::arb,
    x::Int, lim::Float64)
    
    NOverlap = RF(0)
    
    for s in 0 : min(l1,l2)
        NOverlap += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        CTwoCenterOverlap(n1, l1, n2, l2, s, ρ, τ,lim))
    end
    NOverlap
end

function  VTwoCenterOverlap(
    n1::arb, l1::Int, m1::Int,
    n2::arb, l2::Int, m2::Int,
    ρ::arb, τ::arb, theta::arb, phi::arb,
    x::Int, lim::Float64)
        
    NOverlap = RF(0)
        
    for s in 0 : min(l1,l2)
        NOverlap += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        VTwoCenterOverlap(n1, l1, n2, l2, s, ρ, τ,lim))
    end
     NOverlap
end

function  STwoCenterOverlap(
    n1::arb, l1::Int, m1::Int,
    n2::arb, l2::Int, m2::Int,
    ρ::arb, τ::arb, theta::arb, phi::arb,
    x::Int, lim::Float64)
            
    NOverlap = RF(0)
            
    for s in 0 : min(l1,l2)
        NOverlap += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        STwoCenterOverlap(n1, l1, n2, l2, s, ρ, τ,lim))
    end
    NOverlap
end

function  TwoCenterOverlapRec(
    n1::arb, l1::Int, m1::Int,
    n2::arb, l2::Int, m2::Int,
    ρ::arb, τ::arb, theta::arb, phi::arb,
    x::Int, lim::Int)
        NOverlap = RF(0)
    
    for s in 0 : min(l1,l2)
        NOverlap += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        TwoCenterOverlapRec(n1, l1, n2, l2, s, ρ, τ,lim))
    end
        NOverlap
end

function TwoCenterOverlapRec(
    n1::arb, l1::Int, m1:: Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb, x::Int, lim::Int)

    TwoCenterOverlapRec(n1, l1, m1, n2, l2, m2, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), theta, phi, x, lim)
end
###########################################################################################################
###########################################################################################################
################################# TWO CENTER NUCLEAR ATTRACTION INTEGRALS #################################(43)
################################################### ALIGNED MOLECULAR COORDINATE SYSTEM
function  TwoCenterNucAttractABA(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Int)

    coef = (2 * ζ1) // (Sqrt(2 * n1 * (2 * n1 - 1)))
    restco = TwoCenterOverlap(n1 - RF(1), l1, n2, l2, λ, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2),lim)

    return coef * restco
end

function  CTwoCenterNucAttractABA(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64)

    coef = (2 * ζ1) // (Sqrt(2 * n1 * (2 * n1 - 1)))
    restco = CTwoCenterOverlap(n1 - RF(1), l1, n2, l2, λ, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2),lim)

    return coef * restco
end

function  VTwoCenterNucAttractABA(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64)

    coef = (2 * ζ1) // (Sqrt(2 * n1 * (2 * n1 - 1)))
    restco = VTwoCenterOverlap(n1 - RF(1), l1, n2, l2, λ, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2),lim)

    return coef * restco
end

function  STwoCenterNucAttractABA(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64)

    coef = (2 * ζ1) // (Sqrt(2 * n1 * (2 * n1 - 1)))
    restco = STwoCenterOverlap(n1 - RF(1), l1, n2, l2, λ, (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2),lim)

    return coef * restco
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################### NON-ALIGNED MOLECULAR COORDINATE SYSTEM
function  TwoCenterNucAttractABA(
    n1::arb, l1::Int, m1::Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb,
    x::Int, lim::Int)

    coef = (2 * ζ1) // (Sqrt(2 * n1 * (2 * n1 - 1)))
    restco = TwoCenterOverlap(
    n1 - RF(1), l1, m1, 
    n2, l2, m2, 
    (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), theta, phi, x, lim)

    return coef * restco
end

function  CTwoCenterNucAttractABA(
    n1::arb, l1::Int, m1::Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb,
    x::Int, lim::Float64)

    coef = (2 * ζ1) // (Sqrt(2 * n1 * (2 * n1 - 1)))
    restco = CTwoCenterOverlap(
    n1 - RF(1), l1, m1, 
    n2, l2, m2, 
    (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), theta, phi, x, lim)

    return coef * restco
end

function  VTwoCenterNucAttractABA(
    n1::arb, l1::Int, m1::Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb,
    x::Int, lim::Float64)

    coef = (2 * ζ1) // (Sqrt(2 * n1 * (2 * n1 - 1)))
    restco = VTwoCenterOverlap(
    n1 - RF(1), l1, m1, 
    n2, l2, m2, 
    (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), theta, phi, x, lim)

    return coef * restco
end

function  STwoCenterNucAttractABA(
    n1::arb, l1::Int, m1::Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb,
    x::Int, lim::Float64)

    coef = (2 * ζ1) // (Sqrt(2 * n1 * (2 * n1 - 1)))
    restco = STwoCenterOverlap(
    n1 - RF(1), l1, m1, 
    n2, l2, m2, 
    (R//2)*(ζ1+ζ2), (ζ1-ζ2)//(ζ1+ζ2), theta, phi, x, lim)

    return coef * restco
end
###########################################################################################################
function TwoCenterNucAttractAAB(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Int
)
    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)
    n2z = RF(0.000000000000000000000000000000000000000000000001)

    coeff = (R // 2) * SlaterON(n1, n2, ρ, τ)

    sresaab = RF(0)
    for s1 in abs(l1 - l2) : l1 + l2
        for s2 in (EvenQ(l1 + l2 - λ) ? -λ : -λ + 1) : 2 : l1 + l2
            for s3 in 0 : s2
                sresaab1 = Sqrt(2 * s1 + 1)
                sresaab2 = GauntGF(l1, λ, l2, λ, s1)
                sresaab3 = SphPCG(s3, s2, 0, s1, 0, 0, 0, 0)
                sresaab4 = AuxiliaryG(RF(10//10), s3, n1 + n2 - s2 - 1, n2z, RF(10//10), ρ, ρ, lim)
                sresaab += sresaab1 * sresaab2 * sresaab3 * sresaab4
            end
        end
    end

    resaab = coeff * sresaab
    return resaab
end

function CTwoCenterNucAttractAAB(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64
)
    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)
    n2z = RF(0.000000000000000000000000000000000000000000000001)

    coeff = (R // 2) * SlaterON(n1, n2, ρ, τ)

    sresaab = RF(0)
    for s1 in abs(l1 - l2) : l1 + l2
        for s2 in (EvenQ(l1 + l2 - λ) ? -λ : -λ + 1) : 2 : l1 + l2
            for s3 in 0 : s2
                sresaab1 = Sqrt(2 * s1 + 1)
                sresaab2 = GauntGF(l1, λ, l2, λ, s1)
                sresaab3 = SphPCG(s3, s2, 0, s1, 0, 0, 0, 0)
                sresaab4 = CuhreAuxiliaryG(RF(10//10), s3, n1 + n2 - s2 - 1, n2z, RF(10//10), ρ, ρ, lim)
                sresaab += sresaab1 * sresaab2 * sresaab3 * sresaab4
            end
        end
    end

    resaab = coeff * sresaab
    return resaab
end

function VTwoCenterNucAttractAAB(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64
)
    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)
    n2z = RF(0.000000000000000000000000000000000000000000000001)

    coeff = (R // 2) * SlaterON(n1, n2, ρ, τ)

    sresaab = RF(0)
    for s1 in abs(l1 - l2) : l1 + l2
        for s2 in (EvenQ(l1 + l2 - λ) ? -λ : -λ + 1) : 2 : l1 + l2
            for s3 in 0 : s2
                sresaab1 = Sqrt(2 * s1 + 1)
                sresaab2 = GauntGF(l1, λ, l2, λ, s1)
                sresaab3 = SphPCG(s3, s2, 0, s1, 0, 0, 0, 0)
                sresaab4 = VegasAuxiliaryG(RF(10//10), s3, n1 + n2 - s2 - 1, n2z, RF(10//10), ρ, ρ, lim)
                sresaab += sresaab1 * sresaab2 * sresaab3 * sresaab4
            end
        end
    end

    resaab = coeff * sresaab
    return resaab
end

function STwoCenterNucAttractAAB(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64
)
    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)
    n2z = RF(0.000000000000000000000000000000000000000000000001)

    coeff = (R // 2) * SlaterON(n1, n2, ρ, τ)

    sresaab = RF(0)
    for s1 in abs(l1 - l2) : l1 + l2
        for s2 in (EvenQ(l1 + l2 - λ) ? -λ : -λ + 1) : 2 : l1 + l2
            for s3 in 0 : s2
                sresaab1 = Sqrt(2 * s1 + 1)
                sresaab2 = GauntGF(l1, λ, l2, λ, s1)
                sresaab3 = SphPCG(s3, s2, 0, s1, 0, 0, 0, 0)
                sresaab4 = SuaveAuxiliaryG(RF(10//10), s3, n1 + n2 - s2 - 1, n2z, RF(10//10), ρ, ρ, lim)
                sresaab += sresaab1 * sresaab2 * sresaab3 * sresaab4
            end
        end
    end

    resaab = coeff * sresaab
    return resaab
end
###########################################################################################################
#function  TwoCenterNucAttractAAB(
#    n1::arb, l1::Int, m1::Int, ζ1::arb,
#    n2::arb, l2::Int, m2::Int, ζ2::arb,
#    R::arb, theta::arb, phi::arb,
#    x::Int, lim::Int
#    )

#    lowL = abs(l1-l2)
#    upL = l1 + l2
#    n2z = RF(0.000000000000000000000000000000000000000000000001)
#    ρ = (R//2)*(ζ1+ζ2)
#    τ = (ζ1-ζ2)//(ζ1+ζ2)

#    resaab = RF(0)
#    for s1 in lowL : 2 : upL
#        for s2 in -s1 : s1
#            resaab += ChargeDOneC(n1,l1,m1,ζ1,n2,l2,m2,ζ2,n1+n2-1,s1,s2) * 
#            TwoCenterOverlap(n1+n2-1,s1,s2,n2z,0,0,ρ,τ,theta,phi,1,lim)
#        end
#    end
    
#    return resaab
#end

function  TwoCenterNucAttractAAB(
    n1::arb, l1::Int, m1::Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb, lim::Int
    )

    lowL = abs(l1-l2)
    upL = l1 + l2

    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)

    resaab = RF(0)
    for s1 in lowL : 2 : upL
        for s2 in -s1 : s1
            resaab1 = Sqrt((4 * Pi(RF)) // (2 * s1 + 1))
            resaab += (
            GGauntG(l1, m1, l2, m2, s1, s2) * RSphCA(s2, m1, m2) * 
            OneCenterP(s1, n1, n2, ζ1, ζ2, R) * resaab1 * SphericalHarmonicsS(s1, s2, theta, phi)
            )
        end
    end

    return resaab//2
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################# TWO CENTER KINETIC ENERGY INTEGRALS #####################################(43)
################################################### ALIGNED MOLECULAR COORDINATE SYSTEM
function  TwoCenterKinEnergy(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Int
    )

    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)
    
    coeff1 = -Power(ζ2, 2) // 2
    coeff2 = 4 * n2 * Sqrt(Gamma(2*n2 - 1) // Gamma(2*n2 + 1))

    if 2 * n2 - 3 < 0
        coeff3 = 0
        overlap3 = 0
    else
        coeff3 = 4* (n2 + l2) * (n2 - l2 - 1) * Sqrt(Gamma(2*n2 - 3) // Gamma(2*n2 + 1))
        overlap3 = TwoCenterOverlap(n1, l1, ζ1, n2 - 2, l2, ζ2, λ, R, lim)
    end

    overlap1 = TwoCenterOverlap(n1, l1, ζ1, n2, l2, ζ2, λ, R, lim)
    overlap2 = TwoCenterOverlap(n1, l1, ζ1, n2 - 1, l2, ζ2, λ, R, lim)

    res = coeff1 * (
        overlap1 - coeff2 * overlap2 + coeff3 * overlap3
    )

    return res

end

function  CTwoCenterKinEnergy(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64
    )

    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)
    
    coeff1 = -Power(ζ2, 2) // 2
    coeff2 = 4 * n2 * Sqrt(Gamma(2*n2 - 1) // Gamma(2*n2 + 1))

    if 2 * n2 - 3 < 0
        coeff3 = 0
        overlap3 = 0
    else
        coeff3 = 4* (n2 + l2) * (n2 - l2 - 1) * Sqrt(Gamma(2*n2 - 3) // Gamma(2*n2 + 1))
        overlap3 = CTwoCenterOverlap(n1, l1, ζ1, n2 - 2, l2, ζ2, λ, R, lim)
    end

    overlap1 = CTwoCenterOverlap(n1, l1, ζ1, n2, l2, ζ2, λ, R, lim)
    overlap2 = CTwoCenterOverlap(n1, l1, ζ1, n2 - 1, l2, ζ2, λ, R, lim)

    res = coeff1 * (
        overlap1 - coeff2 * overlap2 + coeff3 * overlap3
    )

    return res

end

function  VTwoCenterKinEnergy(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64
    )

    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)
    
    coeff1 = -Power(ζ2, 2) // 2
    coeff2 = 4 * n2 * Sqrt(Gamma(2*n2 - 1) // Gamma(2*n2 + 1))

    if 2 * n2 - 3 < 0
        coeff3 = 0
        overlap3 = 0
    else
        coeff3 = 4* (n2 + l2) * (n2 - l2 - 1) * Sqrt(Gamma(2*n2 - 3) // Gamma(2*n2 + 1))
        overlap3 = VTwoCenterOverlap(n1, l1, ζ1, n2 - 2, l2, ζ2, λ, R, lim)
    end

    overlap1 = VTwoCenterOverlap(n1, l1, ζ1, n2, l2, ζ2, λ, R, lim)
    overlap2 = VTwoCenterOverlap(n1, l1, ζ1, n2 - 1, l2, ζ2, λ, R, lim)

    res = coeff1 * (
        overlap1 - coeff2 * overlap2 + coeff3 * overlap3
    )

    return res

end

function  STwoCenterKinEnergy(
    n1::arb, l1::Int, ζ1::arb,
    n2::arb, l2::Int, ζ2::arb,
    λ::Int, R::arb, lim::Float64
    )

    ρ = (R // 2) * (ζ1 + ζ2)
    τ = (ζ1 - ζ2) // (ζ1 + ζ2)
    
    coeff1 = -Power(ζ2, 2) // 2
    coeff2 = 4 * n2 * Sqrt(Gamma(2*n2 - 1) // Gamma(2*n2 + 1))

    if 2 * n2 - 3 < 0
        coeff3 = 0
        overlap3 = 0
    else
        coeff3 = 4* (n2 + l2) * (n2 - l2 - 1) * Sqrt(Gamma(2*n2 - 3) // Gamma(2*n2 + 1))
        overlap3 = STwoCenterOverlap(n1, l1, ζ1, n2 - 2, l2, ζ2, λ, R, lim)
    end

    overlap1 = STwoCenterOverlap(n1, l1, ζ1, n2, l2, ζ2, λ, R, lim)
    overlap2 = STwoCenterOverlap(n1, l1, ζ1, n2 - 1, l2, ζ2, λ, R, lim)

    res = coeff1 * (
        overlap1 - coeff2 * overlap2 + coeff3 * overlap3
    )

    return res

end
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################### NON-ALIGNED MOLECULAR COORDINATE SYSTEM
function  TwoCenterKinEnergy(
    n1::arb, l1::Int, m1:: Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb, x::Int, lim::Int
    )

    NKinEnergy = RF(0)

    for s in 0 : min(l1,l2)
        NKinEnergy += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        TwoCenterKinEnergy(n1, l1, ζ1, n2, l2, ζ2, s, R, lim))
    end
    NKinEnergy
end

function  CTwoCenterKinEnergy(
    n1::arb, l1::Int, m1:: Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb, x::Int, lim::Float64
    )

    NKinEnergy = RF(0)

    for s in 0 : min(l1,l2)
        NKinEnergy += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        CTwoCenterKinEnergy(n1, l1, ζ1, n2, l2, ζ2, s, R, lim))
    end
    NKinEnergy
end

function  VTwoCenterKinEnergy(
    n1::arb, l1::Int, m1:: Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb, x::Int, lim::Float64
    )

    NKinEnergy = RF(0)

    for s in 0 : min(l1,l2)
        NKinEnergy += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        VTwoCenterKinEnergy(n1, l1, ζ1, n2, l2, ζ2, s, R, lim))
    end
    NKinEnergy
end

function  STwoCenterKinEnergy(
    n1::arb, l1::Int, m1:: Int, ζ1::arb,
    n2::arb, l2::Int, m2::Int, ζ2::arb,
    R::arb, theta::arb, phi::arb, x::Int, lim::Float64
    )

    NKinEnergy = RF(0)

    for s in 0 : min(l1,l2)
        NKinEnergy += (RotaT(s, l1, m1, l2, m2, theta, phi, x) *
        STwoCenterKinEnergy(n1, l1, ζ1, n2, l2, ζ2, s, R, lim))
    end
    NKinEnergy
end
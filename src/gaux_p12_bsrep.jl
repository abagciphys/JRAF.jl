###########################################################################################################
###########################################################################################################
##################### OVERLAP LIKE G VİA SERRIES REPRESENTATION OF BETA FUNCTOINS #########################
###########################################################################################################
###########################################################################################################
export AuxiliaryG
###########################################################################################################
function AuxiliaryG(n1::arb, q::Int, n2::arb, n3::arb, p1::arb, p2::arb, lim::Int)

    two = parent(n1)(2)
    lowl = lowlq = lowl2q = Zero(lim)
    upl = lim
    uplq = q
    upl2q = lim + 2*q

    coef = ((p1^n1) // (Gamma(n1 + One(n1))))

    m1c1 = (n3) // (n2 + 2*q + One(n2))
    m1c2 = (n2 - n3 + 2*q + One(n2) - p2) // (n2 + 2*q + One(n2))

    m2c1 = (n2) // (n3 + 2*q + One(n3))
    m2c2 = (n3 - n2 + 2*q + One(n3) - p2) // (n3 + 2*q + One(n3))

    m1p = zeros(RF, upl2q - lowl2q + 1)
    m2p = zeros(RF, upl2q - lowl2q + 1)
    f1 = zeros(RF, upl2q - lowl2q + 1)
    f2 = zeros(RF, upl2q - lowl2q + 1)

    m11 = zeros(RF, upl2q - lowl2q + 1)
    m12 = zeros(RF, upl2q - lowl2q + 1)
    m21 = zeros(RF, upl2q - lowl2q + 1)
    m22 = zeros(RF, upl2q - lowl2q + 1)

    m1 = zeros(RF, upl2q - lowl2q + 1)
    m2 = zeros(RF, upl2q - lowl2q + 1)

    m1p[1] = Pochhammer(-n2 - 2*q,0)
    f1[1] = RF(1)
    m2p[1] = Pochhammer(-n3,0)
    f2[1] = RF(1)

    m11[1] = real(AuxiliaryGm(n2 + 2*q + two, n3 - One(n3), p2))
    m12[1] = real(AuxiliaryGm(n2 + 2*q + One(n2), n3, p2))
    m21[1] = real(AuxiliaryGm(n3 + 2*q + two, n2 - One(n2), p2))
    m22[1] = real(AuxiliaryGm(n3 + 2*q + One(n3), n2, p2))

    m1[1] = (1//4) * m1c1 * m11[1] + (1//2) * m1c2 * m12[1]
    m2[1] = (1//4) * m2c1 * m21[1] + (1//2) * m2c2 * m22[1]

    for s in lowl2q : upl2q-1
        s1 = s - lowl2q + 1

        m1p[s1+1] = (-n2 - 2*q + s1 - One(n2)) * m1p[s1]
        m2p[s1+1] = (-n3 + s1 -One(n2)) * m2p[s1]
        f1[s1+1] = f2[s1+1] = (s+1) * f1[s1]

        m1c1 = (n3 + s1) // (n2 + 2*q + One(n2) - s1)
        m1c2 = (n2 - n3 - 2*s1 + 2*q + One(n2)-p2) // (n2 + 2*q + One(n2) - s1)
        m2c1 = (n2 + s1) // (n3 + 2*q + One(n3) - s1)
        m2c2 = (n3 - n2 - 2*s1 + 2*q + One(n3)-p2) // (n3 + 2*q + One(n3) - s1)

        m1[s1+1] = (1//4) * m1c1 * m12[s1] + (1//2) * m1c2 * m1[s1]
        m2[s1+1] = (1//4) * m2c1 * m22[s1] + (1//2) * m2c2 * m2[s1]
        m12[s1+1] = m1[s1]
        m22[s1+1] = m2[s1]
    end
    #################
    #################
    #################
    B = zeros(RF, uplq - lowlq + 1)
    b = zeros(RF, uplq - lowlq + 1)

    P1 = zeros(RF, upl - lowl + 1, uplq - lowlq + 1)
    P2 = zeros(RF, upl - lowl + 1, uplq - lowlq + 1)

    M1 = zeros(RF, upl - lowl + 1, uplq - lowlq + 1)
    M2 = zeros(RF, upl - lowl + 1, uplq - lowlq + 1)

    M1P = zeros(RF, upl - lowl + 1, uplq - lowlq + 1)
    M2P = zeros(RF, upl - lowl + 1, uplq - lowlq + 1)

    for s in lowl : upl
        s1 = s - lowl + 1

        P1[s1,1] = m1p[s1] // (n3 + s + 1)
        P2[s1,1] = m2p[s1] // (n2 + 2*q + s + 1)

        M1[s1,1] = m1[s1] // f1[s1]
        M2[s1,uplq+1] = m2[s1] // f2[s1]
    end

    B[1] = Beta(n2 + 2*q + One(n2), n3 + One(n3))
    b[1] = Binomial(q,lowlq)

    for s1 in lowlq : uplq - 1
        s1p = uplq + lowlq - s1 - 1
        s11 = s1 - lowlq + 1
        s11p = s1p - lowlq + 1

        RecB = ((n3+1+2*s11-1) // (n2+2*q+1-2*s11)) * ((n3+1+2*s11-2) // (n2+2*q+1-2*s11+1))
        B[s11 + 1] = RecB * B[s11]
        b[s11 + 1] = ((q - s1) // (s1 + 1)) * b[s11]

        for s2 in 2*s1 : lim + 2*s1
            #s2p = lim + 2*s1 + 2*s1 - s2
            s22 = s2 - 2*s1 + 1
            #s22p = s2p - 2*s1 + 1

            RecP11 = ((-n2 - 2*q + 2*s11 - One(n2) + s22-1)//(-n2 - 2*q + 2*s11 - One(n2)))
            RecP12 = ((-n2 - 2*q + 2*s11 - two + s22-1)//(-n2 - 2*q + 2*s11 - two))
            RecP13 = (n3 + 2*s11 + s22 - 2) // (n3 + 2*s11 + s22)
            P1[s22, s11+1] = RecP11 * RecP12 * RecP13 * P1[s22, s11]

            RecP21 = ((-n3 - 2*s11)//(-n3 - 2*s11 + s22-1))
            RecP22 = ((-n3 - 2*s11 + One(n3))//(-n3 - 2*s11 + One(n3) + s22-1))
            RecP23 = (n2 + 2*q - 2*s11 + s22 + 2) // (n2 + 2*q - 2*s11 + s22)
            P2[s22, s11 + 1] = RecP21 * RecP22 * RecP23 * P2[s22, s11]

            M1[s22, s11 + 1] = m1[s2 + 3] // f1[s22]
            M2[s22, s11p] = m2[s2 + 3] // f2[s22]
        end
    end
    #################
    #################
    #################
    for s1 in lowlq : uplq
        s11 = s1 - lowlq + 1

        b[s11] = ((-1)^s1) * b[s11]
        B[s11] = b[s11] * B[s11]

        for s2 in lowl : upl
            s22 = s2 - lowl + 1

            M1P[s22, s11] =  b[s11] * P1[s22, s11] * M1[s22, s11]
            M2P[s22, s11] =  b[s11] * P2[s22, s11] * M2[s22, s11]
        end
    end

    res1 = Times(
        coef,
        (two^(n2 + n3 + 2*q + One(n2))),
        sum(B),
        (ExpIntegralE(-(n2 + n3 + 2*q + One(n2)), p2))
    )
    res2 = coef * sum(M1P)
    res3 = coef * sum(M2P)
    (res1 - res2 - res3) // (two^(2*q))
end

function AuxiliaryG(n1::Real, q::Int, n2::Real, n3::Real, p1::Real, p2::Real, lim::Int)
    AuxiliaryG(RF(n1), q, RF(n2), RF(n3), RF(p1), RF(p2), lim)
end
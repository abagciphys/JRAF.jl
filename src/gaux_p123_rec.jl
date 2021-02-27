###########################################################################################################
###########################################################################################################
############################## OVERLAP LIKE G USING RECCURRENCE RELATIONS #################################
###########################################################################################################
###########################################################################################################
export AuxiliaryGr
###########################################################################################################
function AuxiliaryGr(n1::arb, q::Int, n2::arb, n3::arb, p1::arb, p2::arb, p3::arb, lim::Int)
    zero = parent(n1)(0)
    one = parent(n1)(1)
    two = parent(n1)(2)
    three = parent(n1)(3)
    four = parent(n1)(4)

    lowq = Zero(lim)
    upq1 = q
    upq2 = 2*q

    Grlrm = zeros(RF, upq2 - lowq + 2, upq2 - lowq + 2)
    Grlrp = zeros(RF, upq2 - lowq + 2, upq2 - lowq + 2)
    Grh = zeros(RF, upq2 - lowq + 2, upq2 - lowq + 2)

    GR = zeros(RF, upq2 - lowq + 1)

    Grlrm[1, 1] = real(AuxiliaryGrlm(n1, RF(n2 + 2*q), n3, p1, p3, p2))
    Grlrm[1, 2] = real(AuxiliaryGrlm(n1, n2 + 2*q, n3 + 10//10, p1, p3, p2))

    Grlrp[1, 1] = real(AuxiliaryGrlp(n1, n3, n2 + 2*q, p1, p3, p2))
    Grlrp[2, 1] = real(AuxiliaryGrlp(n1, n3 + 10//10, n2 + 2*q, p1, p3, p2))

    Grh[1, 1] = real(AuxiliaryGrh(n1, n2 + 2*q, n3, p1, p2, p3, lim))
    Grh[1, 2] = real(AuxiliaryGrh(n1, n2 + 2*q, n3 + 10//10, p1, p2, p3, lim))

    GR[1] = AuxiliaryG(n1, 0, n2 + 2*q, n3, p1, p2, p3, lim)

    #GR[2] = (
    #((p2 + p3) // (p2 - p3)) * ((n3 + 1) // (n2)) * GR[1] +
    #((p2) // (p2 - p3)) * ((1) // (n2)) * (
    #Grlrm[1, 2] - Grlrp[2, 1]
    #) +
    #((p3) // (p2 - p3)) * Grh[1, 2]
    #)

    for s1 in lowq : upq2 - 1
        s2 = s1 - lowq + 1

        Grlrm[s2 + 1, s2 + 1] = (
        ((p2) // (n2 + 2*q - s2 + 1)) * Grlrm[s2, s2 + 1]
        -
        ((n3 + s2) // (n2 + 2*q - s2 + 1)) * Grlrm[s2, s2]
        )

        Grlrm[s2 + 1, s2 + 2] = Grlrm[s2, s2 + 1] - two * Grlrm[s2 + 1, s2 + 1]
        #################
        #################
        Grlrp[s2 + 1, s2 + 1] = (
        ((p2) // (n2 + 2*q - s2 + 1)) * Grlrp[s2 + 1, s2]
        -
        ((n3 + s2) // (n2 + 2*q - s2 + 1)) * Grlrp[s2, s2]
        )

        Grlrp[s2 + 2, s2 + 1] = Grlrp[s2 + 1, s2] + two * Grlrp[s2 + 1, s2 + 1]
        #################
        #################
        Grh[s2 + 1, s2 + 1] = (
        ((p3) // (n2 + 2*q - s2 + 1)) * Grh[s2, s2 + 1]
        +
        ((n3 + s2) // (n2 + 2*q - s2 + 1)) * Grh[s2, s2]
        )

        Grh[s2 + 1, s2 + 2] = -Grh[s2, s2 + 1] + two * Grh[s2 + 1, s2 + 1]
        #################
        #################
        GR[s2 + 1] = (
        ((p2 + p3) // (p2 - p3)) * ((n3 + s2) // (n2 + 2*q - s2 + 1)) * GR[s2] +
        ((p2) // (p2 - p3)) * ((1) // (n2 + 2*q - s2 + 1)) * (
        Grlrm[s2, s2 + 1] - Grlrp[s2 + 1, s2]
        ) +
        ((p3) // (p2 - p3)) * ((1) // (n2 + 2*q - s2 + 1)) * Grh[s2, s2 + 1]
        )
    end

    GRR = zeros(RF, upq1 - lowq + 1)
    for s1 in lowq : upq1
        GRR[s1 + 1] = Power((-10//10),(s1)) * Binomial(q,s1) * GR[2*s1 + 1]
    end
    ((RF(10//10))//(RF(20//10)^(RF(2 * q)))) * sum(GRR)
end
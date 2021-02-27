###########################################################################################################
###########################################################################################################
##################### OVERLAP LIKE G VİA SERRIES REPRESENTATION OF BETA FUNCTOINS #########################
###########################################################################################################
###########################################################################################################
export AuxiliaryG
###########################################################################################################
function AuxiliaryG(n1::arb, q::Int, n2::arb, n3::arb, p1::arb, p2::arb, p3::arb, lim::Int)

    zero = parent(n1)(0)
    one = parent(n1)(1)
    two = parent(n1)(2)
    three = parent(n1)(3)
    four = parent(n1)(4)

    lowq = lowlim = lowlimq = low2limq = low2lim2q = low3lim2q = low3lim3q = low4lim3q = low4lim4q =
    low5lim4q = low5lim5q = 0

    upq = q
    uplim = lim
    uplimq = lim + q
    up2limq = 2*lim + q
    up2lim2q = 2*lim + 2*q
    up3lim2q = 3*lim + 2*q
    up3lim3q = 3*lim + 3*q
    up4lim3q = 4*lim + 3*q
    up4lim4q = 4*lim + 4*q
    up5lim4q = 5*lim + 4*q
    up5lim5q = 5*lim + 5*q
    ###################################################
    ###################################################
    m1c1 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    m1c2 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    m2c1 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    m2c2 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)

    m1r1 = zeros(RF, uplim - lowlim + 1, 2)
    m1r2 = zeros(RF, uplim - lowlim + 1, 2)
    m2r1 = zeros(RF, uplim - lowlim + 1, 2)
    m2r2 = zeros(RF, uplim - lowlim + 1, 2)

    m1 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    m2 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    #################
    #################
    p13rc = zeros(RF, uplim - lowlim + 1, 2)
    p23rc = zeros(RF, uplim - lowlim + 1, 2)

    p13cc = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    p23cc = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)

    p13 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    p23 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    #################
    #################
    s4c1 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    s4c2 = zeros(RF, uplim - lowlim + 1, up5lim4q - low5lim4q + 1)
    #################
    #################
    s4f1 = s4f2 = zeros(RF, up5lim4q - low5lim4q + 1)
    #fs2 is t he the gamma function Gamma(S_2 + 1)#
    fs2 = zeros(RF, uplim - lowlim + 1)
    #################
    #################
    p12 = zeros(RF, up2lim2q - low2lim2q + 1, up3lim2q - low3lim2q + 1)
    p22 = zeros(RF, up2lim2q - low2lim2q + 1, up3lim2q - low3lim2q + 1)
    #################
    #################
    p11 = zeros(RF, uplim - lowlim + 1, uplimq - lowlimq + 1)
    p21 = zeros(RF, uplim - lowlim + 1, uplimq - lowlimq + 1)
    #################
    #################
    b1rc = zeros(RF, uplim - lowlim + 1, 2)
    b2rc = zeros(RF, uplim - lowlim + 1, 2)

    b1cc = zeros(RF, uplim - lowlim + 1, uplimq - lowlimq + 1)
    b2cc = zeros(RF, uplim - lowlim + 1, uplimq - lowlimq + 1)

    B = zeros(RF, uplim - lowlim + 1, uplimq - lowlimq + 1)
    B = zeros(RF, uplim - lowlim + 1, uplimq - lowlimq + 1)
    #################
    #################
    Bins32 = zeros(RF,  uplim - lowlim + 1, uplim - lowlim + 1)
    #################
    #################
    E = zeros(RF, uplim - lowlim + 1)
    ###################################################
    ###################################################
    m1[1,1] = real(AuxiliaryGm(n2 + 2*q, n3 + one, p2))
    m1[2,1] = real(AuxiliaryGm(n2 + 2*q + one, n3 + one, p2))
    m1[1,2] = real(AuxiliaryGm(n2 + 2*q - one, n3 + two, p2))
    m1[2,2] = real(AuxiliaryGm(n2 + 2*q, n3 + two, p2))

    m2[1,1] = real(AuxiliaryGm(n3 + 2*q, n2 + one, p2))
    m2[2,1] = real(AuxiliaryGm(n3 + 2*q + one, n2 + one, p2))
    m2[1,2] = real(AuxiliaryGm(n3 + 2*q - one, n2 + two, p2))
    m2[2,2] = real(AuxiliaryGm(n3 + 2*q, n2 + two, p2))
    #################
    #################
    p13[1,1] = Pochhammer(-(n2 + 2*q), 0)
    p13[2,1] = Pochhammer(-(n2 + 2*q + two), one)
    p13[1,2] = Pochhammer(-(n2 + 2*q) + one, 0)
    p13[2,2] = Pochhammer(-(n2 + 2*q + two) + one, one)

    p23[1,1] = Pochhammer(-(n3 + 2*q), 0)
    p23[2,1] = Pochhammer(-(n3 + 2*q + two), one)
    p23[1,2] = Pochhammer(-(n3 + 2*q) + one, 0)
    p23[2,2] = Pochhammer(-(n3 + 2*q + two) + one, one)
    #################
    #################
    s4c1[1,1] = s4c1[2,1] = (n3 + one)
    s4c1[1,2] = s4c1[2,2] = (n3 + two)

    s4c2[1,1] = s4c2[2,1] = (n2 + one)
    s4c2[1,2] = s4c2[2,2] = (n2 + two)
    #################
    #################
    s4f1[1] = s4f2[1] = fs2[1] = RF(1)
    s4f1[2] = s4f2[2] = fs2[2] = RF(1)
    #################
    #################
    z1 = n2 + 2*q + one
    z2 = n3 + one
    B[1,1] = Beta(z1, z2)
    B[1,2] = Beta(z1 - two, z2 + two)

    B[2,1] = Beta(z1 + two, z2)
    B[2,2] = Beta(z1, z2 + two)
    #################
    #################
    p11[1,1] = Pochhammer(-(n2 + 2*q), 0)
    p11[2,1] = Pochhammer(-(n2 + 2*q + two), one)
    p11[1,2] = Pochhammer(-(n2 + 2*q - two), 0)
    p11[2,2] = Pochhammer(-(n2 + 2*q), one)

    p21[1,1] = Pochhammer(-(n3 + 2*q), 0)
    p21[2,1] = Pochhammer(-(n3 + 2*q + two), one)
    p21[1,2] = Pochhammer(-(n3 + 2*q - two), 0)
    p21[2,2] = Pochhammer(-(n3 + 2*q), one)
    #################
    #################
    Bins32[1,1] = RF(1)
    Bins32[2,1] = RF(1)
    Bins32[1,2] = zero
    Bins32[2,2] = RF(1)
    #################
    #################
    α = (n2 + n3 + 2*q + one)
    E[1] =  ExpIntegralE(-α, p2)
    E[2] =  ExpIntegralE(-(α + one), p2)
    ###################################################
    ###################################################
    for s1 in lowlim : uplim - 2
        s11r = s1 - lowlim + 3

        m1r1[s11r,1] = (-(n2 + s11r - 1) - 2*q + one) // (p2)
        m1r2[s11r,1] = ((n2 + s11r - 1) + n3 + 2*q + p2 + one) // (p2)
        m1r1[s11r,2] = (one - (n2 + s11r - 1) - 2*q + one) // (p2)
        m1r2[s11r,2] = ((n2 + s11r - 1) + n3 + 2*q + p2 + one) // (p2)

        m2r1[s11r,1] = (-(n3 + s11r - 1) - 2*q + one) // (p2)
        m2r2[s11r,1] = ((n3 + s11r - 1) + n2 + 2*q + p2 + one) // (p2)
        m2r1[s11r,2] = (one - (n3 + s11r - 1) - 2*q + one) // (p2)
        m2r2[s11r,2] = ((n3 + s11r - 1) + n2 + 2*q + p2 + one) // (p2)

        m1[s11r,1] = (40//10) * m1r1[s11r,1] * m1[s11r - 2,1] + (20//10) * m1r2[s11r,1] * m1[s11r - 1,1]
        m1[s11r,2] = (40//10) * m1r1[s11r,2] * m1[s11r - 2,2] + (20//10) * m1r2[s11r,2] * m1[s11r - 1,2]

        m2[s11r,1] = (40//10) * m2r1[s11r,1] * m2[s11r - 2,1] + (20//10) * m2r2[s11r,1] * m2[s11r - 1,1]
        m2[s11r,2] = (40//10) * m2r1[s11r,2] * m2[s11r - 2,2] + (20//10) * m2r2[s11r,2] * m2[s11r - 1,2]

        p13rc[s11r,1] = (
            (
                (-(n2 + 2*q) - 2*(s11r - 1))
                //
                (-(n2 + 2*q) - (s11r - 1))
            ) * (-(n2 + 2*q) - 2*(s11r - 1) + 1)
        )
        p23rc[s11r,1] = (
            (
                (-(n3 + 2*q) - 2*(s11r - 1))
                //
                (-(n3 + 2*q) - (s11r - 1))
            ) * (-(n3 + 2*q) - 2*(s11r - 1) + 1)
        )
        p13rc[s11r,2] = (
            (
                (-(n2 - one + 2*q) - 2*(s11r - 1))
                //
                (-(n2 - one + 2*q) - (s11r - 1))
            ) * (-(n2 - one + 2*q) - 2*(s11r - 1) + 1)
        )
        p23rc[s11r,2] = (
            (
                (-(n3 - one + 2*q) - 2*(s11r - 1))
                //
                (-(n3 - one + 2*q) - (s11r - 1))
            ) * (-(n3 - one + 2*q) - 2*(s11r - 1) + 1)
        )

        p13[s11r,1] = p13rc[s11r,1] * p13[s11r - 1,1]
        p13[s11r,2] = p13rc[s11r,2] * p13[s11r - 1,2]

        p23[s11r,1] = p23rc[s11r,1] * p23[s11r - 1,1]
        p23[s11r,2] = p23rc[s11r,2] * p23[s11r - 1,2]

        s4c1[s11r,1] = s4c1[1,1]
        s4c1[s11r,2] = s4c1[1,2]

        s4c2[s11r,1] = s4c2[1,1]
        s4c2[s11r,2] = s4c2[1,2]

        fs2[s11r] = (s11r - 1) * fs2[s11r - 1]

        b1rc[s11r,1] = (
            ((z1 + 2*(s11r - 1) - 1) * (z1 + 2*(s11r - 1) - 2))
            //
            ((z1 + z2 + 2*(s11r - 1) - 1) * (z1 + z2 + 2*(s11r - 1) - 2))
        )

        b2rc[s11r,1] = (
            ((z1 + 2*(s11r - 1) - 3) * (z1 + 2*(s11r - 1) - 4))
            //
            ((z1 + z2 + 2*(s11r - 1) - 1) * (z1 + z2 + 2*(s11r - 1) - 2))
        )

        B[s11r,1] =  b1rc[s11r,1] *  B[s11r - 1,1]
        B[s11r,2] =  b2rc[s11r,1] *  B[s11r - 1,2]

        Bins32[s11r,1] = RF(1)
        Bins32[s11r,2] = RF(s11r - 1)

        E[s11r] = (one // p2) * (Exp(-p2) +  (α + s11r - one) * E[s11r - 1])
    end
    #################
    #################
    for s1 in lowlim : uplim
        s11c = s1 - lowlim + 1

        for s2 in low5lim4q : lim + 4*q + 4*s1 - 2
            s21 = s2 - low5lim4q + 3

            m1c1[s11c,s21] = (n3+(s21-1))//((n2+s1)+2*q+one-(s21-1))
            m1c2[s11c,s21] = ((n2+s1)-n3-2*(s21-1)+2*q+one-p2)//((n2+s1)+2*q+one-(s21-1))

            m2c1[s11c,s21] = (n2+(s21-1))//((n3+s1)+2*q+one-(s21-1))
            m2c2[s11c,s21] = ((n3+s1)-n2-2*(s21-1)+2*q+one-p2)//((n3+s1)+2*q+one-(s21-1))

            m1[s11c,s21] = (10//40)*m1c1[s11c,s21]*m1[s11c,s21-2]+(10//20)*m1c2[s11c,s21]*m1[s11c,s21-1]
            m2[s11c,s21] = (10//40)*m2c1[s11c,s21]*m2[s11c,s21-2]+(10//20)*m2c2[s11c,s21]*m2[s11c,s21-1]

            p13cc[s11c,s21] = ((-(n2 + 2*q) - s1 + (s21 - 1) - one) // (-(n2 + 2*q)-2*s1 + (s21 - 1) - 1))
            p23cc[s11c,s21] = ((-(n3 + 2*q) - s1 + (s21 - 1) - one) // (-(n3 + 2*q)-2*s1 + (s21 - 1) - 1))

            p13[s11c,s21] =  p13cc[s11c,s21] * p13[s11c,s21 - 1]
            p23[s11c,s21] =  p23cc[s11c,s21] * p23[s11c,s21 - 1]

            s4c1[s11c,s21] = (n3 + s2 + two + one)
            s4c2[s11c,s21] = (n2 + s2 + two + one)

            if s21 <= lim + 2*q + 2*s11c - 1
                s4f1[s21] = s4f2[s21] = (s21 - 1) * s4f1[s21 - 1]
            end
            #################
            if s11c > 2 && s21 <= s11c
                Bins32[s11c,s21] = Bins32[s11c - 1,s21 - 1] + Bins32[s11c - 1,s21]
            end
        end
    end
    #################
    #################
    for s1 in low2lim2q : up2lim2q
        s11 = s1 - low2lim2q + 1

        p12[s11,1] = one
        p12[s11,2] = -(n2 + 2*q + lim - s11 + 1)

        p22[s11,1] = one
        p22[s11,2] = -(n3 + 2*q + lim - s11 + 1)

        for s2 in low3lim2q : up3lim2q - 2
            s21 = s2 - low3lim2q + 3

            p12[s11,s21] = (-(n2 + 2*q + lim - s11 + 1) + s21 - 1 - 1) * p12[s11,s21 - 1]

            p22[s11,s21] = (-(n3 + 2*q + lim - s11 + 1) + s21 - 1 - 1) * p22[s11,s21 - 1]

        end
    end
    #################
    #################
    for s1 in lowlim : uplim
        s11 = s1 - lowlim + 1

        p11[s11,1] = p13[s11,1]
        p21[s11,1] = p23[s11,1]

        p11[s11,2] = p13[s11,3]
        p21[s11,2] = p23[s11,3]

        for s2 in lowq : upq + s1 - 2
            s21 = s2 - lowq + 3

            p11[s11,s21] = p13[s11,2*s21-1]
            p21[s11,s21] = p23[s11,2*s21-1]

            B[s11,s21] = (
                ((z2 + 2*s21 - 3) * (z2 + 2*s21 - 4))
                //
                ((z1 + 2*s1 - 2*(s21 - 1)) * (z1 + 2*s1 - 2*(s21 - 1) + 1))
            ) *  B[s11,s21 - 1]
        end
    end
    #################
    #################
    LES1 = zeros(RF, upq - lowq + 1, uplim - lowlim + 1, uplim - lowlim + 1, up3lim2q - low3lim2q + 1)
    LES2 = zeros(RF, upq - lowq + 1, uplim - lowlim + 1, uplim - lowlim + 1)
    BES = zeros(RF, upq - lowq + 1, uplim - lowlim + 1, uplim - lowlim + 1)
    JES1 = zeros(RF, upq - lowq + 1, uplim - lowlim + 1, uplim - lowlim + 1)
    JES2 = zeros(RF, upq - lowq + 1, uplim - lowlim + 1)
    GES1 = zeros(RF, upq - lowq + 1)
    GES2 = RF(0)
    GES3 = RF(0)
    #################
    #################
    for s1 in lowq : upq
        for s2 in lowlim : uplim
            for s3 in lowlim : s2

                BES[s1 + 1,s2 + 1,s3 + 1] =
                ((p3^s2)  // fs2[s2 + 1]) *
                (
                B[s2 + 1,s1 + s3 + 1] *
                E[s2 + 1]
                )

                for s4 in low3lim2q : lim + 2*q + 2*s2

                    LES1[s1 + 1,s2 + 1,s3 + 1,s4 + 1] =
                    ((p3^s2)  // fs2[s2 + 1]) *
                    (
                    (   p11[s2 + 1, s1 + s3 + 1] *
                        m1[s2 + 1,2*s1 + 2*s3 + s4 + 1]
                        //
                        (
                            p13[s2 + 1,2*s1 + 2*s3 + s4 + 1]
                            *
                            s4c1[s2 + 1,2*s1 + 2*s3 + s4 + 1]
                            *
                            s4f1[s4 + 1]
                        ) *
                        p12[lim + 2*s1 - s2 + 2*s3 + 1,s4 + 1]
                    ) + #-
                    (   p21[s2 + 1,(q - s1) + (s2 - s3) + 1]  *
                        m2[s2 + 1,2*(q - s1) + 2*(s2 - s3) + s4 + 1]
                        //
                       (
                            p23[s2 + 1,2*(q - s1) + 2*(s2 - s3) + s4 + 1]
                            *
                            s4c2[s2 + 1,2*(q - s1) + 2*(s2 - s3) + s4 + 1]
                            *
                            s4f2[s4 + 1]
                       ) *
                       p22[lim + 2*(q - s1) - s2 + 2*(s2 - s3) + 1,s4 + 1]
                    )
                    )

                    LES2[s1 + 1,s2 + 1,s3 + 1] += LES1[s1 + 1,s2 + 1,s3 + 1,s4 + 1]

                end

                JES1[s1 + 1,s2 + 1,s3 + 1] =
                Power((20//10), (s2)) * (Power((20//10),(n2 + n3 + 2*q + s2 + 1)) * BES[s1 + 1,s2 + 1,s3 + 1]) -
                Power((20//10),(s2)) * LES2[s1 + 1,s2 + 1,s3 + 1]

                JES2[s1 + 1,s2 + 1] +=
                ((10//10) // Power(20//10, 2*s2)) * Power((-10//10),(s3)) * Bins32[s2 + 1, s3 + 1] *
                JES1[s1 + 1,s2 + 1,s3 + 1]

            end

            GES1[s1 + 1] += Power((-10//10),(s1)) * Bins32[q + 1,s1 + 1] *
            Power(-10//10, s2) * JES2[s1 + 1,s2 + 1]

        end
        GES2 += GES1[s1 + 1]
    end

    GES3 = ((10//10) // ((20//10)^(2*q))) * GES2
    Res = ((p1^(n1))//Gamma(n1+1)) * GES3

    return Res
end
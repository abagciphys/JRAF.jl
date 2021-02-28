###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                        TESTING EXPRESSIONS                                              #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
export NumberQ, IntegerQ, OddQ, EvenQ, RationalQ, RealQ, PositiveQ, NonpositiveQ, NegativeQ, NonnegativeQ,
EqualQ, Rationalize, ARationalize, NO, SignBit, CopySign, FlipSign, Positive, Nonpositive, Negative, 
Nonnegative, Zero, One, Sign, KroneckerDelta, Abs, Min, Max, Plus, Subtract, Times, Divide, Power, Sqrt, 
RSqrt, Surd, Hypot, LogHypot, Exp, Expm1, Ln, Lnp1, Log10, Log, Pi, SqrtPi, LogSqrt2Pi, Sin, SinPi, Cos, 
CosPi, Tan, Cot, Csc, Sec, ArcSin, ArcCos, ArcTan, ArcCot, Sinh, Cosh, Tanh, Coth, ArcSinh, ArcCosh, 
ArcTanh, ArcCoth, Factorial, Gamma, UGamma, LGamma, PGamma, QGamma, Pochhammer, Beta, Binomial, GBinomial, 
ExpIntegralEi, ExpIntegralE, Erf, Erfc, Erfi
###########################################################################################################
################################################## NUMBER #################################################
function NumberQ(nmbrx)
    typeof(nmbrx)<:Number
end
###########################################################################################################
################################################# INTEGER #################################################
function IntegerQ(x::arb)
    return Bool(ccall((:arb_is_int, Nemo.libarb), Cint, (Ref{arb},), x))
 end

IntegerQ(x::Real) = IntegerQ(RF(x))

function IntegerQ(x::acb)
    return Bool(ccall((:acb_is_int, Nemo.libarb), Cint, (Ref{acb},), x))
 end

function IntegerQ(x::Complex)
	xre, xim = reim(x)
	x = CF(xre, xim)
	return IntegerQ(x)
end

IntegerQ(x...) = IntegerQ.((x))
###########################################################################################################
################################################### ODD ###################################################(44)
OddQ(x::arb) = ifelse(IntegerQ(x), rem(convert(Float64, x), 2) != 0, false)
OddQ(x::Real) = ifelse(IntegerQ(x), rem(x, 2) != 0, false)
OddQ(x::acb) = ifelse(IntegerQ(x), rem(real(convert(ComplexF64, x)),2) != 0, false)
OddQ(x::Complex) = ifelse(IntegerQ(x), rem(real(x), 2) != 0, false)
OddQ(x...) = OddQ.((x))
###########################################################################################################
################################################### EVEN ##################################################
EvenQ(x::arb) = ifelse(IntegerQ(x), rem(convert(Float64, x), 2) == 0, false)
EvenQ(x::Real) = ifelse(IntegerQ(x), rem(x, 2) == 0, false)
EvenQ(x::acb) = ifelse(IntegerQ(x), rem(real(convert(ComplexF64, x)),2) == 0, false)
EvenQ(x::Complex) = ifelse(IntegerQ(x), rem(real(x), 2) == 0, false)
EvenQ(x...) = EvenQ.((x))
###########################################################################################################
################################################# RATIONAL ################################################
function RationalQ(rtnlx)
	typeof(rtnlx) == Rational{Int64} ||
	typeof(rtnlx) == Rational{Int128} ||
	typeof(rtnlx) == Rational{BigInt} ||
	typeof(rtnlx) == Complex{Rational{Int64}} ||
	typeof(rtnlx) == Complex{Rational{Int128}} ||
	typeof(rtnlx) == Complex{Rational{BigInt}} ?
	true : false
end
RationalQ(rtnlx...) = RationalQ.((rtnlx))
###########################################################################################################
################################################## REAL ###################################################
function RealQ(x::acb)
    return Bool(ccall((:acb_is_real, Nemo.libarb), Cint, (Ref{acb},), x))
 end

 RealQ(x::Complex) = isreal(x)
###########################################################################################################
################################################# POSITIVE ################################################
################# Return `false` for x<=0 and `true` for x>0.
function PositiveQ(x::arb)
    return Bool(ccall((:arb_is_positive, Nemo.libarb), Cint, (Ref{arb},), x))
 end

PositiveQ(x::Real) = PositiveQ(RF(x))
PositiveQ(x::acb) = RealQ(x) && PositiveQ(real(x))
PositiveQ(x::Complex) = RealQ(x) && PositiveQ(real(x))
###########################################################################################################
############################################### NONPOSITIVE ###############################################
################# Return `true` for x<=0 and `false` for x>0.
function NonpositiveQ(x::arb)
    return Bool(ccall((:arb_is_nonpositive, Nemo.libarb), Cint, (Ref{arb},), x))
 end

NonpositiveQ(x::Real) = NonpositiveQ(RF(x))
NonpositiveQ(x::acb) = RealQ(x) && NonpositiveQ(real(x))
NonpositiveQ(x::Complex) = RealQ(x) && NonpositiveQ(real(x))
###########################################################################################################
################################################# NEGATIVE ################################################
################# Return `false` for x>=0 and `true` for x<0.
function NegativeQ(x::arb)
    return Bool(ccall((:arb_is_negative, Nemo.libarb), Cint, (Ref{arb},), x))
 end

NegativeQ(x::Real) = NegativeQ(RF(x))
NegativeQ(x::acb) = RealQ(x) && NegativeQ(real(x))
NegativeQ(x::Complex) = RealQ(x) && NegativeQ(real(x))
###########################################################################################################
############################################### NONNEGATIVE ###############################################
################# Return `true` for x>=0 and `false` for x<0.
function NonnegativeQ(x::arb)
    return Bool(ccall((:arb_is_nonnegative, Nemo.libarb), Cint, (Ref{arb},), x))
 end

NonnegativeQ(x::Real) = NonnegativeQ(RF(x))
NonnegativeQ(x::acb) = RealQ(x) && NonnegativeQ(real(x))
NonnegativeQ(x::Complex) = RealQ(x) && NonnegativeQ(real(x))
###########################################################################################################
################################################## EQUAL ##################################################
function EqualQ(x::arb, y::arb)
    r = ccall((:arb_equal, Nemo.libarb), Cint, (Ref{arb}, Ref{arb}), x, y)
    return Bool(r)
end
EqualQ(x::Real, y::Real) = EqualQ(RF(x), RF(y))

function EqualQ(x::acb, y::acb)
    r = ccall((:acb_equal, Nemo.libarb), Cint, (Ref{acb}, Ref{acb}), x, y)
    return Bool(r)
end
EqualQ(x::Complex, y::Complex) = EqualQ(CF(real(x), imag(x)), CF(real(y), imag(y)))

EqualQ(x::arb, y::Real) = EqualQ(x, RF(y))
EqualQ(x::Real, y::arb) = EqualQ(RF(x), y)


EqualQ(x::acb, y::Complex) = EqualQ(x, CF(real(y), imag(y)))
EqualQ(x::Complex, y::acb) = EqualQ(CF(real(x), imag(x)), y)

EqualQ(x::arb, y::acb) = ifelse(imag(y) == 0, EqualQ(x, real(y)), false)
EqualQ(x::acb, y::arb) = ifelse(imag(x) == 0, EqualQ(real(x), y), false)

EqualQ(x::arb, y::Complex) = ifelse(imag(y) == 0, EqualQ(x, real(y)), false)
EqualQ(x::Complex, y::arb) = ifelse(imag(x) == 0, EqualQ(real(x), y), false)

EqualQ(x::acb, y::Real) = ifelse(imag(x) == 0, EqualQ(real(x), y), false)
EqualQ(x::Real, y::acb) = ifelse(imag(y) == 0, EqualQ(x, real(y)), false)

EqualQ(x::Real, y::Complex) = ifelse(imag(y) == 0, EqualQ(x, real(y)), false)
EqualQ(x::Complex, y::Real) = ifelse(imag(x) == 0, EqualQ(real(x), y), false)
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                     NUMBER REPRESENTATIONS                                              #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
############################################## ARGUMENT LIST ##############################################
ArgListOpt(Opt,arglox) = arglox

ArgListOpt(Opt,arglox,argloy) = Opt(arglox,argloy)

ArgListOpt(Opt,arglox,argloy,argloz...) = ArgListOpt(Opt, Opt(arglox,argloy), argloz...)
###########################################################################################################
############################################### RATIONALIZE ###############################################(40)
Rationalize(x::Rational) = x
Rationalize(x::Complex{Rational{Int64}}) = x
Rationalize(x::Complex{Rational{Int128}}) = x
Rationalize(x::Complex{Rational{BigInt}}) = x

function Rationalize(x::String)
    if occursin("im", "$x") == false
        return rationalize(BigInt, BigFloat(x))
    elseif occursin("+", "$x") == true
        arrx = split(x, "+")
        restrx = BigFloat(arrx[1])
        if occursin("*", "$arrx[2]") == true
            imstrx = BigFloat(split(arrx[2], "*")[1])
        else
            imstrx = BigFloat(split(arrx[2], "i")[1])
        end
        return rationalize(BigInt, restrx) + im*rationalize(BigInt, imstrx)
    else
        arrx = split(x, "-")
        restrx = BigFloat(arrx[1])
        if occursin("*", "$arrx[2]") == true
            imstrx = BigFloat(split(arrx[2], "*")[1])
        else
            imstrx = BigFloat(split(arrx[2], "i")[1])
        end
        return rationalize(BigInt, restrx) - im*rationalize(BigInt, imstrx)
    end
end

function Rationalize(x::Real)
	restrx = real(x)
	strx = "$restrx"
	return rationalize(BigInt, BigFloat(strx))
end

function Rationalize(x::Complex)
	restrx = BigFloat(real(x))
	imstrx = BigFloat(imag(x))
	reres = rationalize(BigInt, restrx)
	imres = rationalize(BigInt, imstrx)
	return reres + imres*im
end
###########################################################################################################
############################################ ARRAY RATIONALIZE ############################################
ARationalize(x...) = Rationalize.((x))

#function ARationalize(x)
#	arrx = []
#	for s in 1 : length(x)
#		append!(arrx, Rationalize(x[s]))
#	end
#	return arrx
#end
###########################################################################################################
############################################### NEMO OUTPUT ###############################################
function NO(x::arb)
    x = "$x"
    if occursin("+/-", "$x") == true
        x = split("$x", " +/-")[1]
        x = split("$x", "[")[2]
        BigFloat(x)
    else
        BigFloat(x)
    end
end

function NO(x::acb)
    return NO(real(x)) + im * NO(imag(x))
end

NO(x...) = NO.((x))
###########################################################################################################
################################################# SIGNBIT #################################################
################# Returns true if the value of the sign of x is negative, otherwise false.
SignBit(x) = NegativeQ(x)
###########################################################################################################
################################################ COPY SIGN ################################################
################# Return `z` which has the magnitude of `x` and the same sign as `y`.
CopySign(x, y) = ifelse(SignBit(x) != SignBit(y), -x, x)
###########################################################################################################
################################################ FLIP SIGN ################################################
################# Return `x` with its sign flipped if `y` is negative.
FlipSign(x, y) = ifelse(SignBit(y), -x, x)
###########################################################################################################
############################################## POSITIVE SIGN ##############################################(39)
################# Return `-1` for x<=0 and `1` for x>0.
Positive(x) = ifelse(PositiveQ(x), parent(x)(1), parent(x)(-1))
###########################################################################################################
############################################ NONPOSITIVE SIGN #############################################
################# Return `1` for x<=0 and `-1` for x>0.
Nonpositive(x) = ifelse(NonpositiveQ(x), parent(x)(1), parent(x)(-1))
###########################################################################################################
############################################## NEGATIVE SIGN ##############################################
################# Return `-1` for x>=0 and `1` for x<0.
Negative(x) = ifelse(NegativeQ(x), parent(x)(1), parent(x)(-1))
###########################################################################################################
############################################ NONNEGATIVE SIGN #############################################
################# Return `1` for x>=0 and `-1` for x<0.
Nonnegative(x) = ifelse(NonnegativeQ(x), parent(x)(1), parent(x)(-1))
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                          COMPLEX COMPONENTS                                             #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################## UNITY ##################################################
function Zero(r::AcbField)
    z = acb()
    z.parent = r
    return z
end
Zero(R::ArbField) = R(0)

function Zero(z::arb)
    ccall((:arb_zero, Nemo.libarb), Nothing, (Ref{arb},), z)
    return z
end

function Zero(z::acb)
    ccall((:acb_zero, Nemo.libarb), Nothing, (Ref{acb},), z)
    return z
 end
 Zero(x::Number) = oftype(x,0)
 Zero(::Type{T}) where {T<:Number} = convert(T,0)

 function One(r::AcbField)
    z = acb()
    ccall((:acb_one, Nemo.libarb), Nothing, (Ref{acb}, ), z)
    z.parent = r
    return z
  end
One(R::ArbField) = R(1)

One(x::acb) = parent(x)(1)
One(x::arb) = parent(x)(1)
One(::Type{T}) where {T<:Number} = convert(T,1)
One(onex::T) where {T<:Number} = One(T)
###########################################################################################################
################################################### SIGN ##################################################
#Sign(x::Number) = x == 0 ? x / Abs(oneunit(x)) : x / Abs(x)
Sign(x::arb) = ifelse(x < 0, parent(One(x))(-1), ifelse(x > 0, One(x), parent(One(x))(x)))
Sign(x::Real) = ifelse(x < 0, oftype(One(x),-1), ifelse(x > 0, One(x), parent(One(x))(x)))
Sign(x::acb) = x == 0 ? x / Abs(One(x)) : x / Abs(x)
Sign(x::Complex) = x == 0 ? x / Abs(One(x)) : x / Abs(x)

Sign(x, y) = Times(Sign(x), Sign(y))

Sign(a, b, c, xs...) = ArgListOpt(Sign, Sign(Sign(a,b),c), xs...)
###########################################################################################################
############################################# KRONECKER DELTA #############################################(38)
KroneckerDelta(x, y) = ifelse(EqualQ(x,y), One(x), zero(x))
###########################################################################################################
################################################ ABSOLUTE #################################################(41)
function Abs(x::arb)
    z = parent(x)()
    ccall((:arb_abs, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}), z, x)
    return z
end

function Abs(x::Real)
    z = BigFloat()
    ccall((:mpfr_abs, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
    z, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end

function Abs(x::acb)
    z = acb()
    ccall((:acb_abs, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    z.parent = AcbField(parent(x).prec)
    return z
end
Abs(x::Complex) = Abs(CF(real(x), imag(x)))

Abs(x...) = Abs.((x))
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                         ORDER STATISTICS                                                #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################# MINIMUM #################################################
function Min(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_min, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end

function Min(x::Real, y::Real)
	z = BigFloat()
	ccall(
	(:mpfr_min, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, y, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end
Min(x::acb, y::acb) = imag(x) == 0 && imag(y) == 0 && Min(real(x), real(y))
Min(x::Complex, y::Complex) = imag(x) == 0 && imag(y) == 0 && Min(real(x), real(y))

Min(x::arb, y::acb) =  Min(CF(x, parent(x)(0)), y)
Min(x::arb, y::Real) = Min(x, RF(y))
Min(x::arb, y::Complex) = Min(CF(x, parent(x)(0)), CF(real(y), imag(y)))

Min(x::acb, y::arb) =  Min(x, CF(y, parent(y)(0)))
Min(x::acb, y::Real) = Min(x, CF(y, parent(y)(0)))
Min(x::acb, y::Complex) = Min(x, CF(real(y), imag(y)))

Min(x::Real, y::arb) = Min(RF(x), y)
Min(x::Real, y::acb) = Min(CF(x, parent(x)(0)), y)
Min(x::Real, y::Complex) = Min(Complex(x, oftype(x,0)), y)

Min(x::Complex, y::arb) = Min(CF(real(x), imag(x)), CF(y, parent(y)(0)))
Min(x::Complex, y::acb) = Min(CF(real(x), imag(x)), y)
Min(x::Complex, y::Real) = Min(x, Complex(y, oftype(y,0)))

Min(a, b, c, xs...) = ArgListOpt(Min, Min(Min(a,b),c), xs...)
###########################################################################################################
################################################# MAXIMUM #################################################
function Max(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_max, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end

function Max(x::Real, y::Real)
	z = BigFloat()
	ccall(
	(:mpfr_max, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, y, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

Max(x::acb, y::acb) = imag(x) == 0 && imag(y) == 0 && Max(real(x), real(y))
Max(x::Complex, y::Complex) = imag(x) == 0 && imag(y) == 0 && Max(real(x), real(y))

Max(x::arb, y::acb) =  Max(CF(x, parent(x)(0)), y)
Max(x::arb, y::Real) = Max(x, RF(y))
Max(x::arb, y::Complex) = Max(CF(x, parent(x)(0)), CF(real(y), imag(y)))

Max(x::acb, y::arb) =  Max(x, CF(y, parent(y)(0)))
Max(x::acb, y::Real) = Max(x, CF(y, parent(y)(0)))
Max(x::acb, y::Complex) = Max(x, CF(real(y), imag(y)))

Max(x::Real, y::arb) = Max(RF(x), y)
Max(x::Real, y::acb) = Max(CF(x, parent(x)(0)), y)
Max(x::Real, y::Complex) = Max(Complex(x, oftype(x,0)), y)

Max(x::Complex, y::arb) = Max(CF(real(x), imag(x)), CF(y, parent(y)(0)))
Max(x::Complex, y::acb) = Max(CF(real(x), imag(x)), y)
Max(x::Complex, y::Real) = Max(x, Complex(y, oftype(y,0)))

Max(a, b, c, xs...) = ArgListOpt(Max, Max(Max(a,b),c), xs...)
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                         ARITHMETIC FUNCTIONS                                            #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################# ADDITION ################################################(42)
function Plus(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_add, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end
function Plus(x::Real, y::Real)
	z = BigFloat()
	ccall((:mpfr_add, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, y, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Plus(x::acb, y::acb)
    z = parent(x)()
    ccall((:acb_add, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Int),
    z, x, y, parent(x).prec)
    return z
end
Plus(x::Complex, y::Complex) = Plus(CF(real(x), imag(x)), CF(real(y), imag(y)))

Plus(x::arb, y::acb) = Plus(CF(x, parent(x)(0)), y)
Plus(x::arb, y::Real) = Plus(x, RF(y))
Plus(x::arb, y::Complex) = Plus(CF(x, parent(x)(0)), CF(real(y), imag(y)))

Plus(x::acb, y::arb) = Plus(x, CF(y, parent(y)(0)))
Plus(x::acb, y::Real) = Plus(x, CF(y, parent(y)(0)))
Plus(x::acb, y::Complex) = Plus(x, CF(real(y), imag(y)))

Plus(x::Real, y::arb) = Plus(RF(x), y)
Plus(x::Real, y::acb) = Plus(CF(x, parent(x)(0)), y)
Plus(x::Real, y::Complex) = Plus(Complex(x, oftype(x,0)), y)

Plus(x::Complex, y::arb) = Plus(CF(real(x), imag(x)), CF(y, parent(y)(0)))
Plus(x::Complex, y::acb) = Plus(CF(real(x), imag(x)), y)
Plus(x::Complex, y::Real) = Plus(x, Complex(y, oftype(y,0)))

Plus(a, b, c, xs...) = ArgListOpt(Plus, Plus(Plus(a,b),c), xs...)
###########################################################################################################
############################################### SUBTRACTİON ###############################################
function Subtract(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_sub, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end
function Subtract(x::Real, y::Real)
	z = BigFloat()
	ccall((:mpfr_sub, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, y, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Subtract(x::acb, y::acb)
    z = parent(x)()
    ccall((:acb_sub, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Int),
    z, x, y, parent(x).prec)
    return z
end
Subtract(x::Complex, y::Complex) = Subtract(CF(real(x), imag(x)), CF(real(y), imag(y)))

Subtract(x::arb, y::acb) = Subtract(CF(x, parent(x)(0)), y)
Subtract(x::arb, y::Real) = Subtract(x, RF(y))
Subtract(x::arb, y::Complex) = Subtract(CF(x, parent(x)(0)), CF(real(y), imag(y)))

Subtract(x::acb, y::arb) = Subtract(x, CF(y, parent(y)(0)))
Subtract(x::acb, y::Real) = Subtract(x, CF(y, parent(y)(0)))
Subtract(x::acb, y::Complex) = Subtract(x, CF(real(y), imag(y)))

Subtract(x::Real, y::arb) = Subtract(RF(x), y)
Subtract(x::Real, y::acb) = Subtract(CF(x, parent(x)(0)), y)
Subtract(x::Real, y::Complex) = Subtract(Complex(x, oftype(x,0)), y)

Subtract(x::Complex, y::arb) = Subtract(CF(real(x), imag(x)), CF(y, parent(y)(0)))
Subtract(x::Complex, y::acb) = Subtract(CF(real(x), imag(x)), y)
Subtract(x::Complex, y::Real) = Subtract(x, Complex(y, oftype(y,0)))

Subtract(a, b, c, xs...) = ArgListOpt(Subtract, Subtract(Subtract(a,b),c), xs...)
###########################################################################################################
############################################## MULTIPLICATION #############################################
function Times(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_mul, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end
function Times(x::Real, y::Real)
	z = BigFloat()
	ccall((:mpfr_mul, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, y, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Times(x::acb, y::acb)
    z = parent(x)()
    ccall((:acb_mul, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Int),
    z, x, y, parent(x).prec)
    return z
end
Times(x::Complex, y::Complex) = Times(CF(real(x), imag(x)), CF(real(y), imag(y)))

Times(x::arb, y::acb) = Times(CF(x, parent(x)(0)), y)
Times(x::arb, y::Real) = Times(x, RF(y))
Times(x::arb, y::Complex) = Times(CF(x, parent(x)(0)), CF(real(y), imag(y)))

Times(x::acb, y::arb) = Times(x, CF(y, parent(y)(0)))
Times(x::acb, y::Real) = Times(x, CF(y, parent(y)(0)))
Times(x::acb, y::Complex) = Times(x, CF(real(y), imag(y)))

Times(x::Real, y::arb) = Times(RF(x), y)
Times(x::Real, y::acb) = Times(CF(x, parent(x)(0)), y)
Times(x::Real, y::Complex) = Times(Complex(x, oftype(x,0)), y)

Times(x::Complex, y::arb) = Times(CF(real(x), imag(x)), CF(y, parent(y)(0)))
Times(x::Complex, y::acb) = Times(CF(real(x), imag(x)), y)
Times(x::Complex, y::Real) = Times(x, Complex(y, oftype(y,0)))

Times(a, b, c, xs...) = ArgListOpt(Times, Times(Times(a,b),c), xs...)
###########################################################################################################
################################################ DIVISION #################################################
function Divide(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_div, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end
function Divide(x::Real, y::Real)
	z = BigFloat()
	ccall((:mpfr_div, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, y, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Divide(x::acb, y::acb)
    z = parent(x)()
    ccall((:acb_div, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Int),
    z, x, y, parent(x).prec)
    return z
end
Divide(x::Complex, y::Complex) = Divide(CF(real(x), imag(x)), CF(real(y), imag(y)))

Divide(x::arb, y::acb) = Divide(CF(x, parent(x)(0)), y)
Divide(x::arb, y::Real) = Divide(x, RF(y))
Divide(x::arb, y::Complex) = Divide(CF(x, parent(x)(0)), CF(real(y), imag(y)))

Divide(x::acb, y::arb) = Divide(x, CF(y, parent(y)(0)))
Divide(x::acb, y::Real) = Divide(x, CF(y, parent(y)(0)))
Divide(x::acb, y::Complex) = Divide(x, CF(real(y), imag(y)))

Divide(x::Real, y::arb) = Divide(RF(x), y)
Divide(x::Real, y::acb) = Divide(CF(x, parent(x)(0)), y)
Divide(x::Real, y::Complex) = Divide(Complex(x, oftype(x,0)), y)

Divide(x::Complex, y::arb) = Divide(CF(real(x), imag(x)), CF(y, parent(y)(0)))
Divide(x::Complex, y::acb) = Divide(CF(real(x), imag(x)), y)
Divide(x::Complex, y::Real) = Divide(x, Complex(y, oftype(y,0)))

Divide(a, b, c, xs...) = ArgListOpt(Divide, Divide(Divide(a,b),c), xs...)
################# Return the multiplicative inverse of x, i.e. 1/x
function Divide(x::arb)
    z = parent(x)()
    ccall((:arb_inv, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return parent(x)(z)
end
Divide(x::Real) = Divide(RF(x))

function Divide(x::acb)
    z = parent(x)()
    ccall((:acb_inv, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end
Divide(x::Complex) = Divide(CF(real(x), imag(x)))
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
# 									  ELEMENTARY FUNCTIONS                                                #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################## POWER ##################################################(43)
function Power(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_pow, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end
function Power(x::Real, y::Real)
	z = BigFloat()
	ccall((:mpfr_pow, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, y, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Power(x::acb, y::acb)
    z = parent(x)()
    ccall((:acb_pow, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Int),
    z, x, y, parent(x).prec)
    return z
end
Power(x::Complex, y::Complex) = Power(CF(real(x), imag(x)), CF(real(y), imag(y)))

Power(x::arb, y::acb) = Power(CF(x, parent(x)(0)), y)
Power(x::arb, y::Real) = Power(x, RF(y))
Power(x::arb, y::Complex) = Power(CF(x, parent(x)(0)), CF(real(y), imag(y)))

Power(x::acb, y::arb) = Power(x, CF(y, parent(y)(0)))
Power(x::acb, y::Real) = Power(x, CF(y, parent(y)(0)))
Power(x::acb, y::Complex) = Power(x, CF(real(y), imag(y)))

Power(x::Real, y::arb) = Power(RF(x), y)
Power(x::Real, y::acb) = Power(CF(x, parent(x)(0)), y)
Power(x::Real, y::Complex) = Power(Complex(x, oftype(x,0)), y)

Power(x::Complex, y::arb) = Power(CF(real(x), imag(x)), CF(y, parent(y)(0)))
Power(x::Complex, y::acb) = Power(CF(real(x), imag(x)), y)
Power(x::Complex, y::Real) = Power(x, Complex(y, oftype(y,0)))

Power(a, b, c, xs...) = ArgListOpt(Power, Power(Power(a,b),c), xs...)
###########################################################################################################
############################################### SQUARE ROOT ###############################################
function Sqrt(x::arb)
    z = parent(x)()
    ccall((:arb_sqrt, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Sqrt(x::Real)
    z = BigFloat()
    ccall((:mpfr_sqrt, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
    z, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end

function Sqrt(x::acb)
    z = parent(x)()
    ccall((:acb_sqrt, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Sqrt(x::Complex) = Sqrt(CF(real(x), imag(x)))
###########################################################################################################
########################################### INVERSE SQUARE ROOT ###########################################
################# Return 1/Sqrt(x)
function RSqrt(x::arb)
    z = parent(x)()
    ccall((:arb_rsqrt, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function RSqrt(x::Real)
    z = BigFloat()
    ccall((:mpfr_rec_sqrt, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
    z, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end

function RSqrt(x::acb)
    z = parent(x)()
    ccall((:acb_rsqrt, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

RSqrt(x::Complex) = RSqrt(CF(real(x), imag(x)))
###########################################################################################################
################################################ NTH ROOT #################################################(49)
function Surd(x::arb, k::UInt)
    z = parent(x)()
    ccall((:arb_root, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, UInt, Int), z, x, k, parent(x).prec)
    return z
end
Surd(x::arb, k::Int) = Surd(x, UInt(k))

function Surd(x::Real, y::UInt)
    z = BigFloat()
    ccall((:mpfr_root, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, UInt, Int32),
    z, x, y, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end
Surd(x::Real, k::Int) = Surd(x, UInt(k))

function Surd(x::acb, k::UInt)
    z = parent(x)()
    ccall((:acb_root_ui, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, UInt, Int), z, x, k, parent(x).prec)
    return z
end
Surd(x::acb, k::Int) = Surd(x, UInt(k))

Surd(x::Complex, k::UInt) = Surd(CF(real(x), imag(x)), UInt(k))
Surd(x::Complex, k::Int) = Surd(x, UInt(k))
###########################################################################################################
################################################## HYPOT ##################################################
function Hypot(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_hypot, Nemo.libarb), Nothing,(Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end

function Hypot(x::Real, y::Real)
    z = BigFloat()
    ccall((:mpfr_hypot, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
    z, x, y, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end
################# Return Log(Sqrt(x^2 + y^3)).
function LogHypot(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_log_hypot, Nemo.libarb), Nothing,(Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end
LogHypot(x::Real, y::Real) = LogHypot(RF(x), RF(y))
###########################################################################################################
############################################### EXPONENTIAL ###############################################
function Exp(x::arb)
    z = parent(x)()
    ccall((:arb_exp, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Exp(x::Real)
	z = BigFloat()
	ccall((:mpfr_exp, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Exp(x::acb)
    z = parent(x)()
    ccall((:acb_exp, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Exp(x::Complex) = Exp(CF(real(x), imag(x)))
################# Return Exp(x) - 1, evaluated accurately for small x.
function Expm1(x::arb)
    z = parent(x)()
    ccall((:arb_expm1, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
 end

 function Expm1(x::Real)
	z = BigFloat()
	ccall((:mpfr_expm1, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Expm1(x::acb)
    z = parent(x)()
    ccall((:acb_expm1, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
 end
 Expm1(x::Complex) = Expm1(CF(real(x), imag(x)))
###########################################################################################################
################################################ LOGARİTHM ################################################
function Ln(x::arb)
    z = parent(x)()
    ccall((:arb_log, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Ln(x::Real)
	z = BigFloat()
	ccall((:mpfr_log, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Ln(x::acb)
    z = parent(x)()
    ccall((:acb_log, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Ln(x::Complex) = Ln(CF(real(x), imag(x)))
################# Return Ln(1+x), evaluated accurately for small x.
function Lnp1(x::arb)
    z = parent(x)()
    ccall((:arb_log1p, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Lnp1(x::Real)
	z = BigFloat()
	ccall((:mpfr_log1p, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Lnp1(x::acb)
    z = parent(x)()
    ccall((:acb_log1p, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end
Lnp1(x::Complex) = Lnp1(CF(real(x), imag(x)))
################# Log_{10}
function Log10(r::ArbField)
    z = r()
    ccall((:arb_const_log10, Nemo.libarb), Nothing, (Ref{arb}, Int), z, prec(r))
    return z
end

Log(x::arb) = Divide(Ln(x), Log10(ArbField(sprec)))

function Log(x::Real)
	z = BigFloat()
	ccall((:mpfr_log10, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

Log(x::acb) = Divide(Ln(x), Log10(ArbField(sprec)))
Log(x::Complex) = Divide(Ln(x), Log10(ArbField(sprec)))
################# Log_{b}(x), in Mathematica Log(b,x)
function Log(x::arb, b::UInt)
    z = parent(x)()
    ccall((:arb_log_base_ui, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, UInt, Int),
    z, x, b, parent(x).prec)
    return z
end
Log(x::arb, b::Int) = Log(x, UInt(b))
Log(x::arb, b::arb) = Log(x, Int(convert(Float64, b)))
Log(x::Real, b::Int) = Log(RF(x), UInt(b))
Log(x::Real, b::Real) = Log(RF(x), UInt(b))
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                     TRIGONOMETRIC FUNCTIONS                                             #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
function Pi(r::ArbField)
    z = r()
    ccall((:arb_const_pi, Nemo.libarb), Nothing, (Ref{arb}, Int), z, prec(r))
    return z
end

function Pi(r::AcbField)
    z = r()
    ccall((:acb_const_pi, Nemo.libarb), Nothing, (Ref{acb}, Int), z, prec(r))
    return z
end
################# Sqrt(Pi)
function SqrtPi(r::ArbField)
    z = r()
    ccall((:arb_const_sqrt_pi, Nemo.libarb), Nothing, (Ref{arb}, Int), z, prec(r))
    return z
end
################# Log(Sqrt(2Pi))
function LogSqrt2Pi(r::ArbField)
    z = r()
    ccall((:arb_const_log_sqrt2pi, Nemo.libarb), Nothing, (Ref{arb}, Int), z, prec(r))
    return z
end
###########################################################################################################
################################################### SINE ##################################################(44)
function Sin(x::arb)
    z = parent(x)()
    ccall((:arb_sin, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Sin(x::Real)
	z = BigFloat()
	ccall((:mpfr_sin, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Sin(x::acb)
    z = parent(x)()
    ccall((:acb_sin, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Sin(x::Complex) = Sin(CF(real(x), imag(x)))
################# Sin(Pi*x)
function SinPi(x::arb)
    z = parent(x)()
    ccall((:arb_sin_pi, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end
SinPi(x::Real) = SinPi(RF(x))

function SinPi(x::acb)
    z = parent(x)()
    ccall((:acb_sin_pi, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end
SinPi(x::Complex) = SinPi(CF(real(x), imag(x)))
###########################################################################################################
################################################## COSINE #################################################
function Cos(x::arb)
    z = parent(x)()
    ccall((:arb_cos, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
 end

 function Cos(x::Real)
	z = BigFloat()
    ccall((:mpfr_cos, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end

function Cos(x::acb)
    z = parent(x)()
    ccall((:acb_cos, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Cos(x::Complex) = Cos(CF(real(x), imag(x)))
################# Cos(Pi*x)
function CosPi(x::arb)
    z = parent(x)()
    ccall((:arb_cos_pi, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end
CosPi(x::Real) = CosPi(RF(x))

function CosPi(x::acb)
    z = parent(x)()
    ccall((:acb_cos_pi, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end
CosPi(x::Complex) = CosPi(CF(real(x), imag(x)))
###########################################################################################################
################################################# TANGENT #################################################
function Tan(x::arb)
    z = parent(x)()
    ccall((:arb_tan, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Tan(x::Real)
	z = BigFloat()
	ccall((:mpfr_tan, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Tan(x::acb)
    z = parent(x)()
    ccall((:acb_tan, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Tan(x::Complex) = Tan(CF(real(x), imag(x)))
###########################################################################################################
################################################ COTANGENT ################################################
function Cot(x::arb)
    z = parent(x)()
    ccall((:arb_cot, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Cot(x::Real)
	z = BigFloat()
	ccall((:mpfr_cot, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Cot(x::acb)
    z = parent(x)()
    ccall((:acb_cot, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Cot(x::Complex) = Cot(CF(real(x), imag(x)))
###########################################################################################################
################################################# COSECAND ################################################
function Csc(x::arb)
	z = parent(x)()
	ccall((:arb_csc, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
	return z
end

function Csc(x::Real)
	z = BigFloat()
	ccall((:mpfr_csc, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Csc(x::acb)
	z = parent(x)()
	ccall((:acb_csc, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
	return z
end

Csc(x::Complex) = Csc(CF(real(x), imag(x)))
###########################################################################################################
################################################## SECAND #################################################
function Sec(x::arb)
	z = parent(x)()
	ccall((:arb_sec, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
	return z
end

function Sec(x::Real)
	z = BigFloat()
	ccall((:mpfr_cot, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Sec(x::acb)
	z = parent(x)()
	ccall((:acb_sec, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
	return z
end

Sec(x::Complex) = Sec(CF(real(x), imag(x)))
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                             INVERSE TRIGONOMETRIC FUNCTIONS                                             #
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################# ARCSINE #################################################(42)
function ArcSin(x::arb)
	z = parent(x)()
	ccall((:arb_asin, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
	return z
end

function ArcSin(x::Real)
	z = BigFloat()
	ccall((:mpfr_asin, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function ArcSin(x::acb)
	z = parent(x)()
	ccall((:acb_asin, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
	return z
end

ArcSin(x::Complex) = ArcSin(CF(real(x), imag(x)))
###########################################################################################################
################################################ ARCCOSINE ################################################
function ArcCos(x::arb)
	z = parent(x)()
	ccall((:arb_acos, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
	return z
end

function ArcCos(x::Real)
	z = BigFloat()
	ccall((:mpfr_acos, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function ArcCos(x::acb)
	z = parent(x)()
	ccall((:acb_acos, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
	return z
end

ArcCos(x::Complex) = ArcCos(CF(real(x), imag(x)))
###########################################################################################################
################################################ ARCTANGENT ###############################################
function ArcTan(x::arb)
    z = parent(x)()
    ccall((:arb_atan, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function ArcTan(x::Real)
	z = BigFloat()
	ccall((:mpfr_atan, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function ArcTan(x::acb)
    z = parent(x)()
    ccall((:acb_atan, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

ArcTan(x::Complex) = ArcTan(CF(real(x), imag(x)))
################# ArcTan(y,x) = ArcTan(y/x)
function ArcTan(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_atan2, Nemo.libarb), Nothing,
                (Ref{arb}, Ref{arb}, Ref{arb}, Int), z, x, y, parent(x).prec)
    return z
end

function ArcTan(y::Real, x::Real)
    z = BigFloat()
    ccall((:mpfr_atan2, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int),
    z, y, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end
ArcTan(x::acb, y::acb) = imag(x) == 0 && imag(y) == 0 && ArcTan(real(x), real(y))
ArcTan(x::Complex, y::Complex) = imag(x) == 0 && imag(y) == 0 && ArcTan(real(x), real(y))

ArcTan(x::arb, y::acb) =  ArcTan(CF(x, parent(x)(0)), y)
ArcTan(x::arb, y::Real) = ArcTan(x, RF(y))
ArcTan(x::arb, y::Complex) = ArcTan(CF(x, parent(x)(0)), CF(real(y), imag(y)))

ArcTan(x::acb, y::arb) =  ArcTan(x, CF(y, parent(y)(0)))
ArcTan(x::acb, y::Real) = ArcTan(x, CF(y, parent(y)(0)))
ArcTan(x::acb, y::Complex) = ArcTan(x, CF(real(y), imag(y)))

ArcTan(x::Real, y::arb) = ArcTan(RF(x), y)
ArcTan(x::Real, y::acb) = ArcTan(CF(x, parent(x)(0)), y)
ArcTan(x::Real, y::Complex) = ArcTan(Complex(x, oftype(x,0)), y)

ArcTan(x::Complex, y::arb) = ArcTan(CF(real(x), imag(x)), CF(y, parent(y)(0)))
ArcTan(x::Complex, y::acb) = ArcTan(CF(real(x), imag(x)), y)
ArcTan(x::Complex, y::Real) = ArcTan(x, Complex(y, oftype(y,0)))

ArcTan(a, b, c, xs...) = ArgListOpt(ArcTan, ArcTan(ArcTan(a,b),c), xs...)
###########################################################################################################
############################################### ARCCOTANGENT ##############################################
ArcCot(x) = ArcTan(Divide(One(x), x))
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                     HYPERBOLIC FUNCTIONS                                                #
###########################################################################################################
###########################################################################################################
###########################################################################################################
############################################# HYPERBOLIC SINE #############################################(38)
function Sinh(x::arb)
    z = parent(x)()
    ccall((:arb_sinh, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Sinh(x::Real)
	z = BigFloat()
	ccall((:mpfr_sinh, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
    z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Sinh(x::acb)
    z = parent(x)()
    ccall((:acb_sinh, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Sinh(x::Complex) = Sinh(CF(real(x), imag(x)))
###########################################################################################################
############################################ HYPERBOLIC COSINE ############################################
function Cosh(x::arb)
    z = parent(x)()
    ccall((:arb_cosh, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Cosh(x::Real)
	z = BigFloat()
	ccall((:mpfr_cosh, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Cosh(x::acb)
    z = parent(x)()
    ccall((:acb_cosh, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Cosh(x::Complex) = Cosh(CF(real(x), imag(x)))
###########################################################################################################
############################################ HYPERBOLIC TANGENT ###########################################
function Tanh(x::arb)
    z = parent(x)()
    ccall((:arb_tanh, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Tanh(x::Real)
	z = BigFloat()
	ccall((:mpfr_tanh, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Tanh(x::acb)
    z = parent(x)()
    ccall((:acb_tanh, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Tanh(x::Complex) = Tanh(CF(real(x), imag(x)))
###########################################################################################################
########################################### HYPERBOLIC COTANGENT ##########################################
function Coth(x::arb)
    z = parent(x)()
    ccall((:arb_coth, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Coth(x::Real)
	z = BigFloat()
	ccall((:mpfr_coth, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function Coth(x::acb)
    z = parent(x)()
    ccall((:acb_coth, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Coth(x::Complex) = Coth(CF(real(x), imag(x)))
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                  INVERSE HYPERBOLIC FUNCTIONS						      		          #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
############################################ HYPERBOLIC ARCSINE ###########################################(37)
function ArcSinh(x::arb)
	z = parent(x)()
	ccall((:arb_asinh, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
	return z
end

function ArcSinh(x::Real)
	z = BigFloat()
	ccall((:mpfr_asinh, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function ArcSinh(x::acb)
	z = parent(x)()
	ccall((:acb_asinh, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
	return z
end

ArcSinh(x::Complex) = ArcSinh(CF(real(x), imag(x)))
###########################################################################################################
########################################### HYPERBOLIC ARCCOSINE ##########################################
function ArcCosh(x::arb)
	z = parent(x)()
	ccall((:arb_acosh, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
	return z
end

function ArcCosh(x::Real)
	z = BigFloat()
	ccall((:mpfr_acosh, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function ArcCosh(x::acb)
	z = parent(x)()
	ccall((:acb_acosh, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
	return z
end

ArcCosh(x::Complex) = ArcCosh(CF(real(x), imag(x)))
###########################################################################################################
########################################## HYPERBOLIC ARCTANGENT ##########################################
function ArcTanh(x::arb)
	z = parent(x)()
	ccall((:arb_atanh, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
	return z
end

function ArcTanh(x::Real)
	z = BigFloat()
	ccall(
	(:mpfr_atanh, :libmpfr),Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function ArcTanh(x::acb)
	z = parent(x)()
	ccall((:acb_atanh, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
	return z
end

ArcTanh(x::Complex) = ArcTanh(CF(real(x), imag(x)))
###########################################################################################################
######################################### HYPERBOLIC ARCCOTANGENT #########################################
ArcCoth(x) = ArcTanh(Divide(One(x), x))
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#                              			GAMMA BETA ERF								      	              #
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################################ FACTORİALS ###############################################
function Factorial(n::UInt, r::ArbField)
    z = r()
    ccall((:arb_fac_ui, Nemo.libarb), Nothing, (Ref{arb}, UInt, Int), z, n, r.prec)
    return z
end
Factorial(n::Int, r::ArbField) = n < 0 ? Factorial(r(n)) : Factorial(UInt(n), r)

Factorial(x::arb) = Gamma(Plus(x, One(x)))
Factorial(x::Real) = Factorial(RF(x))
Factorial(x::acb) = Gamma(Plus(x, One(x)))
Factorial(x::Complex) = Factorial(CF(real(x), imag(x)))
###########################################################################################################
################################################## GAMMA ##################################################
function Gamma(x::arb)
    z = parent(x)()
    ccall((:arb_gamma, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Gamma(x::Real)
    z = BigFloat()
    ccall((:mpfr_gamma, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end

function Gamma(x::acb)
    z = parent(x)()
    ccall((:acb_gamma, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

Gamma(x::Complex) = Gamma(CF(real(x), imag(x)))
################# Log(Gamma(x))

###########################################################################################################
######################################### UPPER INCOMPLETE GAMMA ##########################################
function UGamma(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_hypgeom_gamma_upper, Nemo.libarb), Nothing,(Ref{arb}, Ref{arb}, Ref{arb}, Int, Int),
    z, x, y, 0, parent(x).prec)
    return z
end

function UGamma(x::Real, y::Real)
	z = BigFloat()
	ccall((:mpfr_gamma_inc, :libmpfr),Int32,(Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, y, Base.MPFR.ROUNDING_MODE[])
	return RF(z)
end

function UGamma(x::acb, y::acb)
    z = parent(x)()
    ccall((:acb_hypgeom_gamma_upper, Nemo.libarb), Nothing,(Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
    z, x, y, 0, parent(x).prec)
    return z
end
UGamma(x::Complex, y::Complex) = UGamma(CF(real(x), imag(x)), CF(real(x), imag(x)))

UGamma(x::arb, y::acb) =  UGamma(CF(x, parent(x)(0)), y)
UGamma(x::arb, y::Real) = UGamma(x,RF(y))
UGamma(x::arb, y::Complex) = UGamma(CF(x, parent(x)(0)), CF(real(y), imag(y)))

UGamma(x::acb, y::arb) =  UGamma(x, CF(y, parent(y)(0)))
UGamma(x::acb, y::Real) = UGamma(x, CF(y, parent(y)(0)))
UGamma(x::acb, y::Complex) = UGamma(x, CF(real(y), imag(y)))

UGamma(x::Real, y::arb) = UGamma(RF(x),y)
UGamma(x::Real, y::acb) = UGamma(CF(x, parent(x)(0)), y)
UGamma(x::Real, y::Complex) = UGamma(CF(x, parent(x)(0)), CF(real(y), imag(y)))

UGamma(x::Complex, y::arb) = UGamma(CF(real(x), imag(x)), CF(y, parent(y)(0)))
UGamma(x::Complex, y::acb) = UGamma(CF(real(x), imag(x)), y)
UGamma(x::Complex, y::Real) = UGamma(CF(real(x), imag(x)), CF(y, parent(y)(0)))

UGamma(a, b, c, xs...) = ArgListOpt(UGamma, UGamma(UGamma(a,b),c), xs...)
###########################################################################################################
######################################### LOWER INCOMPLETE GAMMA ##########################################
function LGamma(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_hypgeom_gamma_lower, Nemo.libarb), Nothing,(Ref{arb}, Ref{arb}, Ref{arb}, Int, Int),
    z, x, y, 0, parent(x).prec)
    return z
end
LGamma(x::Real, y::Real) = LGamma(RF(x), RF(y))

function LGamma(x::acb, y::acb)
    z = parent(x)()
    ccall((:acb_hypgeom_gamma_lower, Nemo.libarb), Nothing,(Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
    z, x, y, 0, parent(x).prec)
    return z
end
LGamma(x::Complex, y::Complex) = LGamma(CF(real(x), imag(x)), CF(real(x), imag(x)))

LGamma(x::arb, y::acb) =  LGamma(CF(x, parent(x)(0)), y)
LGamma(x::arb, y::Real) = LGamma(x,RF(y))
LGamma(x::arb, y::Complex) = LGamma(CF(x, parent(x)(0)), CF(real(y), imag(y)))

LGamma(x::acb, y::arb) =  LGamma(x, CF(y, parent(y)(0)))
LGamma(x::acb, y::Real) = LGamma(x, CF(y, parent(y)(0)))
LGamma(x::acb, y::Complex) = LGamma(x, CF(real(y), imag(y)))

LGamma(x::Real, y::arb) = LGamma(RF(x),y)
LGamma(x::Real, y::acb) = LGamma(CF(x, parent(x)(0)), y)
LGamma(x::Real, y::Complex) = LGamma(CF(x, parent(x)(0)), CF(real(y), imag(y)))

LGamma(x::Complex, y::arb) = LGamma(CF(real(x), imag(x)), CF(y, parent(y)(0)))
LGamma(x::Complex, y::acb) = LGamma(CF(real(x), imag(x)), y)
LGamma(x::Complex, y::Real) = LGamma(CF(real(x), imag(x)), CF(y, parent(y)(0)))

LGamma(a, b, c, xs...) = ArgListOpt(LGamma, LGamma(LGamma(a,b),c), xs...)
###########################################################################################################
############################################ REGULARIZED GAMMA ############################################
PGamma(x, y) = Divide(LGamma(x, y), Gamma(x))
PGamma(a, b, c, xs...) = ArgListOpt(PGamma, PGamma(PGamma(a,b),c), xs...)
#################
QGamma(x, y) = Divide(UGamma(x, y), Gamma(x))
QGamma(a, b, c, xs...) = ArgListOpt(QGamma, QGamma(QGamma(a,b),c), xs...)
###########################################################################################################
################################################ POCHHAMMER ###############################################
function Pochhammer(x::arb, n::UInt)
    z = parent(x)()
    ccall((:arb_rising_ui, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, UInt, Int),
    z, x, n, parent(x).prec)
    return z
end

Pochhammer(x::arb, n::Int) = NonnegativeQ(n) && Pochhammer(x, UInt(n))
Pochhammer(x::arb, n::arb) = Pochhammer(x, Int(convert(Float64, n)))
Pochhammer(x::Real, n::Int) = NonnegativeQ(n) && Pochhammer(RF(x), UInt(n))
Pochhammer(x::Real, n::Real) = Pochhammer(RF(x), RF(n))

function Pochhammer(x::acb, n::UInt)
    z = parent(x)()
    ccall((:acb_rising_ui, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, UInt, Int),
    z, x, n, parent(x).prec)
    return z
end

Pochhammer(x::acb, n::Int) = NonnegativeQ(n) && Pochhammer(x, UInt(n))
Pochhammer(x::acb, n::acb) = Pochhammer(x, Int(convert(ComplexF64, n)))
Pochhammer(x::Complex, n::Int) = NonnegativeQ(n) && Pochhammer(CF(real(x), imag(x)), UInt(n))
Pochhammer(x::Complex, n::Complex) = Pochhammer(CF(real(x), imag(x)), CF(real(n), imag(n)))

Pochhammer(x::arb, y::acb) =  Pochhammer(CF(x, parent(x)(0)), y)
Pochhammer(x::arb, y::Real) = Pochhammer(x, RF(y))
Pochhammer(x::arb, y::Complex) = Pochhammer(CF(x, parent(x)(0)), CF(real(y), imag(y)))

Pochhammer(x::acb, y::arb) =  Pochhammer(x, CF(y, parent(y)(0)))
Pochhammer(x::acb, y::Real) = Pochhammer(x, CF(y, parent(y)(0)))
Pochhammer(x::acb, y::Complex) = Pochhammer(x, CF(real(y), imag(y)))

Pochhammer(x::Real, y::arb) = Pochhammer(RF(x), y)
Pochhammer(x::Real, y::acb) = Pochhammer(CF(x, parent(x)(0)), y)
Pochhammer(x::Real, y::Complex) = Pochhammer(Complex(x, oftype(x,0)), y)

Pochhammer(x::Complex, y::arb) = Pochhammer(CF(real(x), imag(x)), CF(y, parent(y)(0)))
Pochhammer(x::Complex, y::acb) = Pochhammer(CF(real(x), imag(x)), y)
Pochhammer(x::Complex, y::Real) = Pochhammer(x, Complex(y, oftype(y,0)))
###########################################################################################################
################################################## BETA ###################################################(44)
################# 0 < x <= 1 ! while x = 1, Beta(a,b,x) = Beta(a,b) in mathematica Beta(x,a,b)
function Beta(a::arb, b::arb, x::arb)
    z = parent(x)()
    ccall((:arb_hypgeom_beta_lower, Nemo.libarb), Nothing,(Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Int, Int),
    z, a, b, x, 0, parent(x).prec)
    return z
end
Beta(a::Real, b::Real, x::Real) = Beta(RF(a), RF(b), RF(x))

function Beta(a::acb, b::acb, x::acb)
    z = parent(x)()
    ccall((:arb_hypgeom_beta_lower, Nemo.libarb), Nothing,(Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
    z, a, b, x, 0, parent(x).prec)
    return z
end
Beta(a::Complex, b::Complex, x::Complex) = Beta(
    CF(real(a), imag(a)), CF(real(b), imag(b)), CF(real(x), imag(x))
)

Beta(a, b) = Beta(a, b, parent(a)(1))
###########################################################################################################
################################################ BINOMIAL #################################################
function Binomial(x::arb, n::UInt)
    z = parent(x)()
    ccall((:arb_bin_ui, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, UInt, Int),
    z, x, n, parent(x).prec)
    return z
end

Binomial(x::arb, n::Int) = n < 0 ? 0 : Binomial(x, UInt(n))
Binomial(x::arb, n::arb) = n < 0 ? 0 : Binomial(x, Int(convert(Float64, n)))
Binomial(x::Real, n::Int) = n < 0 ? 0 : Binomial(RF(x), n)
Binomial(x::Real, n::Real) = n < 0 ? 0 : Binomial(RF(x), Int(n))
###########################################################################################################
########################################## GENERALİZED BINOMIAL ###########################################
################# F[n,n',m,N]
function GBinomial(z::Int, x::Int, y::Int)
    llim = Int(((z - x) + abs(z - x)) / 2)
    ulim = min(z, y)
    r = zeros(RF, ulim - llim + 1)
    @threads for s in llim : ulim
        r[s - llim + 1] = Times(Power(-1, s), Binomial(x, z - s), Binomial(y, s))
    end
    return sum(r)
end
GBinomial(x::Real, y::Real, z::Real) = GBinomial(Int(x), Int(y), Int(z))
GBinomial(x::arb, y::arb, z::arb) = GBinomial(NO(x), NO(y), NO(z))
function GBinomial(x::Real, y::Real, z::Int, lim::Int)
		r = zeros(RF, lim + 1)
        @threads for s in 0 : lim
            if s <= z
			    sign = Power(-1, s)
                r[s + 1] = Times(sign, Binomial(x, z - s), Binomial(y, s))
            else
                break
            end
        end
        return sum(r)
end
GBinomial(x::Real, y::Real, z::Real, lim::Real) = GBinomial(x, y, Int(z), Int(lim))
GBinomial(x::arb, y::arb, z::arb, lim::arb) = GBinomial(NO(x), NO(y), NO(z), NO(lim))
###########################################################################################################
########################################## EXPONENTIAL INTEGRALS ##########################################
################# E1(x) = -Ei(-x)
function ExpIntegralEi(x::arb)
    z = parent(x)()
    ccall((:arb_hypgeom_ei, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function ExpIntegralEi(x::Real)
    z = BigFloat()
    ccall((:mpfr_eint, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end

function ExpIntegralEi(x::acb)
    z = parent(x)()
    ccall((:arb_hypgeom_ei, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end

ExpIntegralEi(x::Complex) = ExpIntegralEi(CF(real(x), imag(x)))
################# E_{x}(y)
function ExpIntegralE(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_hypgeom_expint, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Ref{arb}, Int),
    z, x, y, parent(x).prec)
    return z
end
ExpIntegralE(x::Real, y::Real) = ExpIntegralE(RF(x), RF(y))

function ExpIntegralE(x::acb, y::acb)
    z = parent(x)()
    ccall((:acb_hypgeom_expint, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Ref{acb}, Int),
    z, x, y, parent(x).prec)
    return z
end
ExpIntegralE(x::Complex, y::Complex) = ExpIntegralE(CF(real(x), imag(x)), CF(real(y), imag(y)))

ExpIntegralE(x::arb, y::acb) =  ExpIntegralE(CF(x, parent(x)(0)), y)
ExpIntegralE(x::arb, y::Real) = ExpIntegralE(x, RF(y))
ExpIntegralE(x::arb, y::Complex) = ExpIntegralE(CF(x, parent(x)(0)), CF(real(y), imag(y)))

ExpIntegralE(x::acb, y::arb) =  ExpIntegralE(x, CF(y, parent(y)(0)))
ExpIntegralE(x::acb, y::Real) = ExpIntegralE(x, CF(y, parent(y)(0)))
ExpIntegralE(x::acb, y::Complex) = ExpIntegralE(x, CF(real(y), imag(y)))

ExpIntegralE(x::Real, y::arb) = ExpIntegralE(RF(x), y)
ExpIntegralE(x::Real, y::acb) = ExpIntegralE(CF(x, parent(x)(0)), y)
ExpIntegralE(x::Real, y::Complex) = ExpIntegralE(Complex(x, oftype(x,0)), y)

ExpIntegralE(x::Complex, y::arb) = ExpIntegralE(CF(real(x), imag(x)), CF(y, parent(y)(0)))
ExpIntegralE(x::Complex, y::acb) = ExpIntegralE(CF(real(x), imag(x)), y)
ExpIntegralE(x::Complex, y::Real) = ExpIntegralE(x, Complex(y, oftype(y,0)))
###########################################################################################################
########################################## PROBABILITY INTEGRALS ##########################################
function Erf(x::arb)
    z = parent(x)()
    ccall((:arb_hypgeom_erf, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Erf(x::Real)
    z = BigFloat()
    ccall((:mpfr_erf, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end

function Erf(x::acb)
    z = parent(x)()
    ccall((:acb_hypgeom_erf, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end
Erf(x::Complex) = Erf(CF(real(x), imag(x)))
################# Erfc(x) gives the complementary error function 1 - Erf(x).
################# This function avoids catastrophic cancellation for large positive x.
function Erfc(x::arb)
    z = parent(x)()
    ccall((:arb_hypgeom_erfc, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end

function Erfc(x::Real)
    z = BigFloat()
    ccall((:mpfr_erfc, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
	z, x, Base.MPFR.ROUNDING_MODE[])
    return RF(z)
end

function Erfc(x::acb)
    z = parent(x)()
    ccall((:acb_hypgeom_erfc, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end
Erfc(x::Complex) = Erfc(CF(real(x), imag(x)))
################# Erfi(x) gives the imaginary error function -iErf(ix).
function Erfi(x::arb)
    z = parent(x)()
    ccall((:arb_hypgeom_erfi, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}, Int), z, x, parent(x).prec)
    return z
end
Erfi(x::Real) = Erfi(RF(x))

function Erfi(x::acb)
    z = parent(x)()
    ccall((:acb_hypgeom_erfi, Nemo.libarb), Nothing, (Ref{acb}, Ref{acb}, Int), z, x, parent(x).prec)
    return z
end
Erfi(x::Complex) = Erfi(CF(real(x), imag(x)))
###########################################################################################################
################################################### END ###################################################
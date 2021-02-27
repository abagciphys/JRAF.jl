using JRAF
using Test
###########################################################################################################
###########################################################################################################
###########################################################################################################
##################################### SIMPLE TEST FOR DATA TYPES ##########################################
x170 = Divide(RF(2), RF(3));
x171 = Divide(3, 2);
x172 = Divide(CF(2,0), CF(3,0));
x173 = Divide(3+0im,2+0im);
x174 = Divide(RF(2), CF(3,0));
x175 = Divide(RF(3), 2);
x176 = Divide(RF(2), 3+0im);
x177 = Divide(CF(3,0), RF(2));
x178 = Divide(CF(2,0), 3);
x179 = Divide(CF(3,0), 2+0im);
x180 = Divide(2, RF(3));
x181 = Divide(3, CF(2,0));
x182 = Divide(2, 3+0im);
x183 = Divide(2+0im, RF(3));
x184 = Divide(3+0im, CF(2,0));
x185 = Divide(CF(3,0), 2);
x186 = Divide(CF(3,0), 2, RF(4), 5+0im);
x187 = Divide(RF(2));
x188 = Divide(2);
x189 = Divide(CF(2,0));
x190 = Divide(2+0im);
###########################################################################################################
x191 = Power(RF(2), RF(3));
x192 = Power(3, 2);
x193 = Power(CF(2,0), CF(3,0));
x194 = Power(3+0im,2+0im);
x195 = Power(RF(2), CF(3,0));
x196 = Power(RF(3), 2);
x197 = Power(RF(2), 3+0im);
x198 = Power(CF(3,0), RF(2));
x199 = Power(CF(2,0), 3);
x200 = Power(CF(3,0), 2+0im);
x201 = Power(2, RF(3));
x202 = Power(3, CF(2,0));
x203 = Power(2, 3+0im);
x204 = Power(2+0im, RF(3));
x205 = Power(3+0im, CF(2,0));
x206 = Power(CF(3,0), 2);
x207 = Power(CF(3,0), 2, RF(4), 5+0im);
###########################################################################################################
x208 = Sqrt(RF(3.5))
x209 = Sqrt(3.5)
x210 = Sqrt(CF(RF(2.3),RF(3)))
x211 = Sqrt(3.5+3im)
###########################################################################################################
x212 = RSqrt(RF(3.5))
x213 = RSqrt(3.5)
x214 = RSqrt(CF(RF(2.3),RF(3)))
x215 = RSqrt(3.5+3im)
###########################################################################################################
###########################################################################################################
@testset "simple test for data types" begin
	@test IntegerQ(2)==true														# 1
	@test IntegerQ(2.2)==false													# 2
	@test IntegerQ(RF(2))==true													# 3
	@test IntegerQ(RF(2.5))==false												# 4
	@test IntegerQ(CF(2,0))==true												# 5
	@test IntegerQ(CF(2,1))==false												# 6
	@test IntegerQ(2+0im)==true													# 7
	@test IntegerQ(2+2im)==false												# 8
	@test IntegerQ(2+2im,2,2.2,RF(2),CF(2,2))==(false,true,false,true,false)	# 9
    @test RealQ(CF(2,0))==true													# 10
	@test RealQ(CF(2,1))==false													# 11
	@test RealQ(2+0im)==true													# 12
	@test RealQ(2+2im)==false													# 13
    @test PositiveQ(CF(2,0))==true												# 14
	@test PositiveQ(CF(-2,0))==false											# 15
	@test PositiveQ(2)==true													# 16
	@test PositiveQ(-2)==false													# 17
	@test PositiveQ(RF(2))==true												# 18
	@test PositiveQ(RF(-2))==false												# 19
	@test PositiveQ(RF(0))==false												# 20
	@test PositiveQ(0)==false													# 21
    @test NonpositiveQ(CF(2,0))==false											# 22
	@test NonpositiveQ(CF(-2,0))==true											# 23
	@test NonpositiveQ(2)==false												# 24
	@test NonpositiveQ(-2)==true												# 25
	@test NonpositiveQ(RF(2))==false											# 26
	@test NonpositiveQ(RF(-2))==true											# 27
	@test NonpositiveQ(RF(0))==true												# 28
	@test NonpositiveQ(0)==true													# 29
    @test NegativeQ(CF(2,0))==false												# 30
	@test NegativeQ(CF(-2,0))==true												# 31
	@test NegativeQ(2)==false													# 32
	@test NegativeQ(-2)==true													# 33
	@test NegativeQ(RF(2))==false												# 34
	@test NegativeQ(RF(-2))==true												# 35
	@test NegativeQ(RF(0))==false												# 36
	@test NegativeQ(0)==false													# 37
    @test NonnegativeQ(CF(2,0))==true											# 38
	@test NonnegativeQ(CF(-2,0))==false											# 39
	@test NonnegativeQ(2)==true													# 40
	@test NonnegativeQ(-2)==false												# 41
	@test NonnegativeQ(RF(2))==true												# 42
	@test NonnegativeQ(RF(-2))==false											# 43
	@test NonnegativeQ(RF(0))==true												# 44
	@test NonnegativeQ(0)==true													# 45
    @test EqualQ(RF(2),RF(2))==true												# 46
	@test EqualQ(RF(2),RF(3))==false											# 47
	@test EqualQ(2, 2)==true													# 48
	@test EqualQ(2, 3)==false													# 49
	@test EqualQ(CF(2), CF(2))==true											# 50
	@test EqualQ(CF(2), CF(3))==false											# 51
	@test EqualQ(2+0im, 2+0im)==true											# 52
	@test EqualQ(2+0im, 3+1im)==false											# 53
	@test EqualQ(RF(2), 2)==true												# 54
	@test EqualQ(RF(2), 3)==false												# 55
	@test EqualQ(2, RF(2))==true												# 56
	@test EqualQ(2, RF(3))==false												# 57
	@test EqualQ(CF(2), 2+0im)==true											# 58
	@test EqualQ(CF(2), 2+1im)==false											# 59
	@test EqualQ(2+0im, CF(2))==true											# 60
	@test EqualQ(2+0im, CF(2,-1))==false										# 61
	@test EqualQ(RF(2), CF(2))==true											# 62
	@test EqualQ(RF(2), CF(2,1))==false											# 63
	@test EqualQ(CF(2), RF(2))==true											# 64
	@test EqualQ(CF(2), RF(3))==false											# 65
	@test EqualQ(RF(2), 2+0im)==true											# 66
	@test EqualQ(RF(2), 2+1im)==false											# 67
	@test EqualQ(2+0im, RF(2))==true											# 68
	@test EqualQ(2+0im, RF(3))==false											# 69
	@test EqualQ(CF(2), 2)==true												# 70
	@test EqualQ(CF(2), 3)==false												# 71
	@test EqualQ(2, CF(2))==true												# 72
	@test EqualQ(2, CF(2,2))==false												# 73
	@test EqualQ(2, 2+0im)==true												# 74
	@test EqualQ(2, 2+1im)==false												# 75
	@test EqualQ(2+0im, 2)==true												# 76
	@test EqualQ(2+0im, 3)==false												# 77
    @test Zero(RF)==RF(0)														# 78
	@test Zero(CF)==CF(0)														# 79
    @test Abs(RF(-5))==5														# 80
	@test Abs(-5)==5															# 81
	@test Abs(CF(-5,0))==5														# 82
	@test Abs(-5+0im)==5														# 83
	@test Abs(RF(-5),5,CF(-5,0),-5+0im)==(5,5,5,5)								# 84
    @test Min(RF(2), RF(3))==2													# 85
	@test Min(3, 2)==2															# 86
	@test Min(CF(2,0), CF(3,0))==2												# 87
	@test Min(3+0im,2+0im)==2													# 88
	@test Min(RF(2), CF(3,0))==2												# 89
	@test Min(RF(3), 2)==2														# 90
	@test Min(RF(2), 3+0im)==2													# 91
	@test Min(CF(3,0), RF(2))==2												# 92
	@test Min(CF(2,0), 3)==2													# 93
	@test Min(CF(3,0), 2+0im)==2												# 94
	@test Min(2, RF(3))==2														# 95
	@test Min(3, CF(2,0))==2													# 96
	@test Min(2, 3+0im)==2														# 97
	@test Min(2+0im, RF(3))==2													# 98
	@test Min(3+0im, CF(2,0))==2												# 99
	@test Min(CF(3,0), 2)==2													# 100
	@test Min(CF(3,0), 2, RF(4), 5+0im)==2										# 101
    @test Max(RF(2), RF(3))==3													# 102
	@test Max(3, 2)==3															# 103
	@test Max(CF(2,0), CF(3,0))==3												# 104
	@test Max(3+0im,2+0im)==3													# 105
	@test Max(RF(2), CF(3,0))==3												# 106
	@test Max(RF(3), 2)==3														# 107
	@test Max(RF(2), 3+0im)==3													# 108
	@test Max(CF(3,0), RF(2))==3												# 109
	@test Max(CF(2,0), 3)==3													# 110
	@test Max(CF(3,0), 2+0im)==3												# 111
	@test Max(2, RF(3))==3														# 112
	@test Max(3, CF(2,0))==3													# 113
	@test Max(2, 3+0im)==3														# 114
	@test Max(2+0im, RF(3))==3													# 115
	@test Max(3+0im, CF(2,0))==3												# 116
	@test Max(CF(3,0), 2)==3													# 117
	@test Max(CF(3,0), 2, RF(4), 5+0im)==5										# 118
    @test Plus(RF(2), RF(3))==5													# 119
	@test Plus(3, 2)==5															# 120
	@test Plus(CF(2,0), CF(3,0))==5												# 121
	@test Plus(3+0im,2+0im)==5													# 122
	@test Plus(RF(2), CF(3,0))==5												# 123
	@test Plus(RF(3), 2)==5														# 124
	@test Plus(RF(2), 3+0im)==5													# 125
	@test Plus(CF(3,0), RF(2))==5												# 126
	@test Plus(CF(2,0), 3)==5													# 127
	@test Plus(CF(3,0), 2+0im)==5												# 128
	@test Plus(2, RF(3))==5														# 129
	@test Plus(3, CF(2,0))==5													# 130
	@test Plus(2, 3+0im)==5														# 131
	@test Plus(2+0im, RF(3))==5													# 132
	@test Plus(3+0im, CF(2,0))==5												# 133
	@test Plus(CF(3,0), 2)==5													# 134
	@test Plus(CF(3,0), 2, RF(4), 5+0im)==14									# 135
    @test Subtract(RF(2), RF(3))==-1											# 136
	@test Subtract(3, 2)==1														# 137
	@test Subtract(CF(2,0), CF(3,0))==-1										# 138
	@test Subtract(3+0im,2+0im)==1												# 139
	@test Subtract(RF(2), CF(3,0))==-1											# 140
	@test Subtract(RF(3), 2)==1													# 141
	@test Subtract(RF(2), 3+0im)==-1											# 142
	@test Subtract(CF(3,0), RF(2))==1											# 143
	@test Subtract(CF(2,0), 3)==-1												# 144
	@test Subtract(CF(3,0), 2+0im)==1											# 145
	@test Subtract(2, RF(3))==-1												# 146
	@test Subtract(3, CF(2,0))==1												# 147
	@test Subtract(2, 3+0im)==-1												# 148
	@test Subtract(2+0im, RF(3))==-1											# 149
	@test Subtract(3+0im, CF(2,0))==1											# 150
	@test Subtract(CF(3,0), 2)==1												# 151
	@test Subtract(CF(3,0), 2, RF(4), 5+0im)==-8								# 152
    @test Times(RF(2), RF(3))==6												# 153
	@test Times(3, 2)==6														# 154
	@test Times(CF(2,0), CF(3,0))==6											# 155
	@test Times(3+0im,2+0im)==6													# 156
	@test Times(RF(2), CF(3,0))==6												# 157
	@test Times(RF(3), 2)==6													# 158
	@test Times(RF(2), 3+0im)==6												# 159
	@test Times(CF(3,0), RF(2))==6												# 160
	@test Times(CF(2,0), 3)==6													# 161
	@test Times(CF(3,0), 2+0im)==6												# 162
	@test Times(2, RF(3))==6													# 163
	@test Times(3, CF(2,0))==6													# 164
	@test Times(2, 3+0im)==6													# 165
	@test Times(2+0im, RF(3))==6												# 166
	@test Times(3+0im, CF(2,0))==6												# 167
	@test Times(CF(3,0), 2)==6													# 168
	@test Times(CF(3,0), 2, RF(4), 5+0im)==120									# 169
    @test NO(Divide(RF(2), RF(3))) == NO(x170)									# 170
	@test NO(Divide(3, 2)) == NO(x171)											# 171
	@test NO(Divide(CF(2,0), CF(3,0))) == NO(x172)								# 172
	@test NO(Divide(3+0im,2+0im)) == NO(x173)									# 173
	@test NO(Divide(RF(2), CF(3,0))) == NO(x174)								# 174
	@test NO(Divide(RF(3), 2)) == NO(x175)										# 175
	@test NO(Divide(RF(2), 3+0im)) == NO(x176)									# 176
	@test NO(Divide(CF(3,0), RF(2))) == NO(x177)								# 177
	@test NO(Divide(CF(2,0), 3)) == NO(x178)									# 178
	@test NO(Divide(CF(3,0), 2+0im)) == NO(x179)								# 179
	@test NO(Divide(2, RF(3))) == NO(x180)										# 180
	@test NO(Divide(3, CF(2,0))) == NO(x181)									# 181
	@test NO(Divide(2, 3+0im)) == NO(x182)										# 182
	@test NO(Divide(2+0im, RF(3))) == NO(x183)									# 183
	@test NO(Divide(3+0im, CF(2,0))) == NO(x184)								# 184
	@test NO(Divide(CF(3,0), 2)) == NO(x185)									# 185
	@test NO(Divide(CF(3,0), 2, RF(4), 5+0im)) == NO(x186)						# 186
	@test NO(Divide(RF(2))) == NO(x187)											# 187
	@test NO(Divide(2)) == NO(x188)												# 188
	@test NO(Divide(CF(2,0))) == NO(x189)										# 189
	@test NO(Divide(2+0im)) == NO(x190)											# 190

    @test NO(Power(RF(2), RF(3))) == NO(x191)									# 191
	@test NO(Power(3, 2)) == NO(x192)											# 192
	@test NO(Power(CF(2,0), CF(3,0))) == NO(x193)								# 193
	@test NO(Power(3+0im,2+0im)) == NO(x194)									# 194
	@test NO(Power(RF(2), CF(3,0))) == NO(x195)									# 195
	@test NO(Power(RF(3), 2)) == NO(x196)										# 196
	@test NO(Power(RF(2), 3+0im)) == NO(x197)									# 197
	@test NO(Power(CF(3,0), RF(2))) == NO(x198)									# 198
	@test NO(Power(CF(2,0), 3)) == NO(x199)										# 199
	@test NO(Power(CF(3,0), 2+0im)) == NO(x200)									# 200
	@test NO(Power(2, RF(3))) == NO(x201)										# 201
	@test NO(Power(3, CF(2,0))) == NO(x202)										# 202
	@test NO(Power(2, 3+0im)) == NO(x203)										# 203
	@test NO(Power(2+0im, RF(3))) == NO(x204)									# 204
	@test NO(Power(3+0im, CF(2,0))) == NO(x205)									# 205
	@test NO(Power(CF(3,0), 2)) == NO(x206)										# 206
	@test NO(Power(CF(3,0), 2, RF(4), 5+0im)) == NO(x207)						# 207

    @test NO(Sqrt(RF(3.5))) == NO(x208)											# 208
	@test NO(Sqrt(3.5)) == NO(x209)												# 209
	@test NO(Sqrt(CF(RF(2.3),RF(3)))) == NO(x210)								# 210
	@test NO(Sqrt(3.5+3im)) == NO(x211)											# 211

    @test NO(RSqrt(RF(3.5))) == NO(x212)										# 212
	@test NO(RSqrt(3.5)) == NO(x213)											# 213
	@test NO(RSqrt(CF(RF(2.3),RF(3)))) == NO(x214)								# 214
	@test NO(RSqrt(3.5+3im)) == NO(x215)										# 215

end
###########################################################################################################
###########################################################################################################
###########################################################################################################
########################################### TEST MATH.JL ##################################################
x216 = Factorial(RF(3.5));
x217 = Factorial(3.5);
x218 = Factorial(CF(RF(2.3),RF(3)));
x219 = Factorial(3.5+3im);
###########################################################################################################
x220 = Gamma(RF(3.5));
x221 = Gamma(3.5);
x222 = Gamma(CF(RF(2.3),RF(3)));
x223 = Gamma(3.5+3im);

###########################################################################################################
x224 = UGamma(RF(2), RF(3));
x225 = UGamma(3, 2);
x226 = UGamma(CF(2,0), CF(3,0));
x227 = UGamma(3+0im,2+0im);
x228 = UGamma(RF(2), CF(3,0));
x229 = UGamma(RF(3), 2);
x230 = UGamma(RF(2), 3+0im);
x231 = UGamma(CF(3,0), RF(2));
x232 = UGamma(CF(2,0), 3);
x233 = UGamma(CF(3,0), 2+0im);
x234 = UGamma(2, RF(3));
x235 = UGamma(3, CF(2,0));
x236 = UGamma(2, 3+0im);
x237 = UGamma(2+0im, RF(3));
x238 = UGamma(3+0im, CF(2,0));
x239 = UGamma(CF(3,0), 2);
x240 = UGamma(CF(3,0), 2, RF(4), 5+0im);

###########################################################################################################
x241 = LGamma(RF(2), RF(3));
x242 = LGamma(3, 2);
x243 = LGamma(CF(2,0), CF(3,0));
x244 = LGamma(3+0im,2+0im);
x245 = LGamma(RF(2), CF(3,0));
x246 = LGamma(RF(3), 2);
x247 = LGamma(RF(2), 3+0im);
x248 = LGamma(CF(3,0), RF(2));
x249 = LGamma(CF(2,0), 3);
x250 = LGamma(CF(3,0), 2+0im);
x251 = LGamma(2, RF(3));
x252 = LGamma(3, CF(2,0));
x253 = LGamma(2, 3+0im);
x254 = LGamma(2+0im, RF(3));
x255 = LGamma(3+0im, CF(2,0));
x256 = LGamma(CF(3,0), 2);
x257 = LGamma(CF(3,0), 2, RF(4), 5+0im);

###########################################################################################################
x258 = Pochhammer(RF(2), RF(3));
x259 = Pochhammer(3, 2);
x260 = Pochhammer(CF(2,0), CF(3,0));
x261 = Pochhammer(3+0im,2+0im);
x262 = Pochhammer(RF(2), CF(3,0));
x263 = Pochhammer(RF(3), 2);
x264 = Pochhammer(RF(2), 3+0im);
x265 = Pochhammer(CF(3,0), RF(2));
x266 = Pochhammer(CF(2,0), 3);
x267 = Pochhammer(CF(3,0), 2+0im);
x268 = Pochhammer(2, RF(3));
x269 = Pochhammer(3, CF(2,0));
x270 = Pochhammer(2, 3+0im);
x271 = Pochhammer(2+0im, RF(3));
x272 = Pochhammer(3+0im, CF(2,0));
x273 = Pochhammer(CF(3,0), 2);

###########################################################################################################
x274 = Beta(RF(3), RF(2), RF(2));
x275 = Beta(3, 2, 2);
x276 = Beta(CF(3,0), CF(2,0), CF(2,0));

###########################################################################################################
x277 = Binomial(RF(3), 2);
x278 = Binomial(RF(3), RF(2));
x279 = Binomial(3.0, 2);
x280 = Binomial(3.0, 2.0);

###########################################################################################################
x281 = GBinomial(3, 2, 3);
x282 = GBinomial(3.0, 2.0, 3.0);
x283 = GBinomial(RF(3), RF(2), RF(3));

###########################################################################################################
x284 = ExpIntegralEi(RF(3));
x285 = ExpIntegralEi(3.0);
x286 = ExpIntegralEi(CF(3,0));
x287 = ExpIntegralEi(3.0 + 0im);

###########################################################################################################
x288 = ExpIntegralE(RF(2), RF(3));
x289 = ExpIntegralE(3, 2);
x290 = ExpIntegralE(CF(2,0), CF(3,0));
x291 = ExpIntegralE(3+0im,2+0im);
x292 = ExpIntegralE(RF(2), CF(3,0));
x293 = ExpIntegralE(RF(3), 2);
x294 = ExpIntegralE(RF(2), 3+0im);
x295 = ExpIntegralE(CF(3,0), RF(2));
x296 = ExpIntegralE(CF(2,0), 3);
x297 = ExpIntegralE(CF(3,0), 2+0im);
x298 = ExpIntegralE(2, RF(3));
x299 = ExpIntegralE(3, CF(2,0));
x300 = ExpIntegralE(2, 3+0im);
x301 = ExpIntegralE(2+0im, RF(3));
x302 = ExpIntegralE(3+0im, CF(2,0));
x303 = ExpIntegralE(CF(3,0), 2);

###########################################################################################################
###########################################################################################################
@testset "math.jl" begin

    @test NO(Factorial(RF(3.5))) == NO(x216)									# 216
    @test NO(Factorial(3.5)) == NO(x217)										# 217
    @test NO(Factorial(CF(RF(2.3),RF(3)))) == NO(x218)							# 218
    @test NO(Factorial(3.5+3im)) == NO(x219)                                    # 219
	@test NO(Factorial(6)) == NO((1//(6+1))*Factorial(7))						# 

    @test NO(Gamma(RF(3.5))) == NO(x220)									    # 220
    @test NO(Gamma(3.5)) == NO(x221)										    # 221
    @test NO(Gamma(CF(RF(2.3),RF(3)))) == NO(x222)							    # 222
    @test NO(Gamma(3.5+3im)) == NO(x223)                                        # 223
	@test Float64(NO(Gamma(3.5))) == Float64(NO((Gamma(3.5 + 1)) // (3.5)))		#

    @test NO(UGamma(RF(2), RF(3))) == NO(x224)									# 224
	@test NO(UGamma(3, 2)) == NO(x225)											# 225
	@test NO(UGamma(CF(2,0), CF(3,0))) == NO(x226)								# 226
	@test NO(UGamma(3+0im,2+0im)) == NO(x227)									# 227
	@test NO(UGamma(RF(2), CF(3,0))) == NO(x228)								# 228
	@test NO(UGamma(RF(3), 2)) == NO(x229)										# 229
	@test NO(UGamma(RF(2), 3+0im)) == NO(x230)									# 230
	@test NO(UGamma(CF(3,0), RF(2))) == NO(x231)								# 231
	@test NO(UGamma(CF(2,0), 3)) == NO(x232)									# 232
	@test NO(UGamma(CF(3,0), 2+0im)) == NO(x233)								# 233
	@test NO(UGamma(2, RF(3))) == NO(x234)										# 234
	@test NO(UGamma(3, CF(2,0))) == NO(x235)									# 235
	@test NO(UGamma(2, 3+0im)) == NO(x236)										# 236
	@test NO(UGamma(2+0im, RF(3))) == NO(x237)									# 237
	@test NO(UGamma(3+0im, CF(2,0))) == NO(x238)								# 238
	@test NO(UGamma(CF(3,0), 2)) == NO(x239)									# 239
	@test NO(UGamma(CF(3,0), 2, RF(4), 5+0im)) == NO(x240)						# 240
    @test NO(LGamma(RF(2), RF(3))) == NO(x241)									# 241
	@test NO(LGamma(3, 2)) == NO(x242)											# 242
	@test NO(LGamma(CF(2,0), CF(3,0))) == NO(x243)								# 243
	@test NO(LGamma(3+0im,2+0im)) == NO(x244)									# 244
	@test NO(LGamma(RF(2), CF(3,0))) == NO(x245)								# 245
	@test NO(LGamma(RF(3), 2)) == NO(x246)										# 246
	@test NO(LGamma(RF(2), 3+0im)) == NO(x247)									# 247
	@test NO(LGamma(CF(3,0), RF(2))) == NO(x248)								# 248
	@test NO(LGamma(CF(2,0), 3)) == NO(x249)									# 249
	@test NO(LGamma(CF(3,0), 2+0im)) == NO(x250)								# 250
	@test NO(LGamma(2, RF(3))) == NO(x251)										# 251
	@test NO(LGamma(3, CF(2,0))) == NO(x252)									# 252
	@test NO(LGamma(2, 3+0im)) == NO(x253)										# 253
	@test NO(LGamma(2+0im, RF(3))) == NO(x254)									# 254
	@test NO(LGamma(3+0im, CF(2,0))) == NO(x255)								# 255
	@test NO(LGamma(CF(3,0), 2)) == NO(x256)									# 256
	@test NO(LGamma(CF(3,0), 2, RF(4), 5+0im)) == NO(x257)						# 257
	@test Float64(NO(UGamma(3,2) + LGamma(3,2))) == Float64(NO(Gamma(3)))		#

    @test NO(Pochhammer(RF(2), RF(3))) == NO(x258)								# 258
	@test NO(Pochhammer(3, 2)) == NO(x259)										# 259
	@test NO(Pochhammer(CF(2,0), CF(3,0))) == NO(x260)							# 260
	@test NO(Pochhammer(3+0im,2+0im)) == NO(x261)								# 261
	@test NO(Pochhammer(RF(2), CF(3,0))) == NO(x262)							# 262
	@test NO(Pochhammer(RF(3), 2)) == NO(x263)									# 263
	@test NO(Pochhammer(RF(2), 3+0im)) == NO(x264)								# 264
	@test NO(Pochhammer(CF(3,0), RF(2))) == NO(x265)							# 265
	@test NO(Pochhammer(CF(2,0), 3)) == NO(x266)								# 266
	@test NO(Pochhammer(CF(3,0), 2+0im)) == NO(x267)							# 267
	@test NO(Pochhammer(2, RF(3))) == NO(x268)									# 268
	@test NO(Pochhammer(3, CF(2,0))) == NO(x269)								# 269
	@test NO(Pochhammer(2, 3+0im)) == NO(x270)									# 270
	@test NO(Pochhammer(2+0im, RF(3))) == NO(x271)								# 271
	@test NO(Pochhammer(3+0im, CF(2,0))) == NO(x272)							# 272
	@test NO(Pochhammer(CF(3,0), 2)) == NO(x273)								# 273
	@test Float64(NO(Pochhammer(3,2))) == Float64(NO(Gamma(3+2) // Gamma(3)))	#

	@test NO(Beta(RF(3), RF(2), RF(2))) == NO(x274)								# 274
	@test NO(Beta(3, 2, 2)) == NO(x275)											# 275
	@test NO(Beta(CF(3,0), CF(2,0), CF(2,0))) == NO(x276)						# 276
	@test Float64(NO(Beta(3, 2, 2))) == Float64(								#
		NO(																		#
			Beta(4, 2, 2) + Beta(3, 3, 2)										#
		)																		#
	)																			#

	@test NO(Binomial(RF(3), 2)) == NO(x277)									# 277
	@test NO(Binomial(RF(3), RF(2))) == NO(x278)								# 278
	@test NO(Binomial(3.0, 2)) == NO(x279)										# 279
	@test NO(Binomial(3.0, 2.0)) == NO(x280)									# 280
	@test Float64(NO(Binomial(7,2))) == Float64(								#
		NO(																		#
			Pochhammer(6,2)// Factorial(2)										#
			)																	#
		)																		#

	@test NO(GBinomial(3, 2, 3)) == NO(x281)									# 281
	@test NO(GBinomial(3.0, 2.0, 3.0)) == NO(x282)								# 282
	@test NO(GBinomial(RF(3), RF(2), RF(3))) == NO(x283)						# 283
	@test Float64(NO(GBinomial(5,10,6))) == Float64(36.0)

	@test NO(ExpIntegralEi(RF(3))) == NO(x284)									# 284
	@test NO(ExpIntegralEi(3.0)) == NO(x285)									# 285
	@test NO(ExpIntegralEi(CF(3,0))) == NO(x286)								# 286
	@test NO(ExpIntegralEi(3.0 + 0im)) == NO(x287)								# 287
	@test Float64(NO(ExpIntegralEi(1))) == Float64(1.8951178163559368)			#

	@test NO(ExpIntegralE(RF(2), RF(3))) == NO(x288)							# 288
	@test NO(ExpIntegralE(3, 2)) == NO(x289)									# 289
	@test NO(ExpIntegralE(CF(2,0), CF(3,0))) == NO(x290)						# 290
	@test NO(ExpIntegralE(3+0im,2+0im)) == NO(x291)								# 291
	@test NO(ExpIntegralE(RF(2), CF(3,0))) == NO(x292)							# 292
	@test NO(ExpIntegralE(RF(3), 2)) == NO(x293)								# 293
	@test NO(ExpIntegralE(RF(2), 3+0im)) == NO(x294)							# 294
	@test NO(ExpIntegralE(CF(3,0), RF(2))) == NO(x295)							# 295
	@test NO(ExpIntegralE(CF(2,0), 3)) == NO(x296)								# 296
	@test NO(ExpIntegralE(CF(3,0), 2+0im)) == NO(x297)							# 297
	@test NO(ExpIntegralE(2, RF(3))) == NO(x298)								# 298
	@test NO(ExpIntegralE(3, CF(2,0))) == NO(x299)								# 299
	@test NO(ExpIntegralE(2, 3+0im)) == NO(x300)								# 300
	@test NO(ExpIntegralE(2+0im, RF(3))) == NO(x301)							# 301
	@test NO(ExpIntegralE(3+0im, CF(2,0))) == NO(x302)							# 302
	@test NO(ExpIntegralE(CF(3,0), 2)) == NO(x303)								# 303
	@test Float64(NO(3 * ExpIntegralE(3+1,2))) == Float64(						#
		NO(																		#
			Exp(-2) - 2 * ExpIntegralE(3,2)										#
		)																		#
		)																		#
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
################################## TEST ANGULAR_COEFFICIENTS.JL ###########################################
smalljlistg = 0:1:2
largejlistg = 0:1:1000
smalljlist = 0:1//2:10
largejlist = 0:1//2:1000

@time @testset "Guseinov's formulae for Clebsch-Gordan coefficients" begin
    for j1 in smalljlistg, j2 in smalljlistg
        d1::Int = 2*j1+1
        d2::Int = 2*j2+1
        M = zeros(Float64, (d1*d2, d1*d2))
        ind1 = 1
        for m1 in -j1:j1, m2 in -j2:j2
            ind2 = 1
            @inbounds for j3 in abs(j1-j2):(j1+j2), m3 in -j3:j3
                M[ind1,ind2] = ClebschGordanG(j1,m1,j2,m2,j3,m3)
                ind2 += 1
            end
            ind1 += 1
        end
        @test M'*M ≈ one(M)
    end
end

@time @testset "Guseinov's formulae for Gaunt coefficients" begin
	for j1 in smalljlistg, j2 in smalljlistg, j3 in abs(j1-j2):(j1+j2)
		for m1 in -j1:j1, m2 in -j2:j2
			g1 = Float64(
				NO(
					GauntG(j1, m1, j2, m2, j3))
					)

			g2 = Float64(
						NO(
							GauntG(j2, m2, j1, m1, j3)
							)
						)

			isnan(g1) ? g1l = 0 : g1l = g1 
			isnan(g2) ? g2l = 0 : g2l = g2

			@test g1l == g2l
		end
	end
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
##################################### TEST SPECİAL_FUNCTIONS.JL ###########################################
lpvar = 0.5
@time @testset "Associated Legendre poloynomials" begin
	for l in smalljlistg, m in l
		lp1 = Float64(NO(LegendreP(l,m,lpvar,0)))
		lp2 = Float64(Plm(l, m, lpvar))
		@test lp1 == lp2
	end
end

@time @testset "Normalized associated Legendre poloynomials" begin
	for l in smalljlistg, m in 0:l
		lp1 = Float32(NO(NLegendreP(l,m,lpvar,0)))
		lp2 = Float32(NO(Nlm(l, m) * Plm(l, m, lpvar) * Sqrt(2*Pi(RF))))
		lp3 = Float32(NO(LegendrePG(l,m,CF(lpvar))))
		@test lp1 == lp2
		@test lp1 == lp3
	end
end
###########################################################################################################
###########################################################################################################
Y1 = computeYlm(pi/2, 0, lmax = smalljlistg[length(smalljlistg)]);
Y2 = computeYlm(
	pi/3, pi/3, 
	lmax = smalljlistg[length(smalljlistg)], 
	SHType = SphericalHarmonics.RealHarmonics()
	);

#SphericalHarmonics.sphericalharmonic(
#	π/3, π/3, l = 500, m = 251,
#	SHType = SphericalHarmonics.RealHarmonics()
#	)

@time @testset "Complex spherical harmonics" begin
	for l in 0:1:2, m in l
		sph1 = Real(Float32(NO(SphericalHarmonicsY(l,m,Pi(RF)//2,RF(0)))))
		sph2 = Real(Float32(Y1[(l,m)]))
		sph3 = Real(Float32(NO(((-1)^m) * SphericalHarmonicsYG(l,m,Pi(RF)//2,RF(0)))))
		@test sph1 == sph2
		@test sph1 == sph3
	end
end

@time @testset "Real spherical harmonics" begin
	for l in 0:1:2, m in l
		sph1 = ((-1)^m) * Float32(NO(real(SphericalHarmonicsSG(l,m,Pi(RF)//3,Pi(RF)//3)+CF(0,0))))
		sph2 = Float32(Y2[(l,m)])
		sph3 = Float32(NO(real(SphericalHarmonicsS(l,m,Pi(RF)//3,Pi(RF)//3)+CF(0,0))))
		@test sph1 == sph2
		@test sph1 == sph3
	end
end

@time @testset "Spherical Spinors" begin
	β = -1
	for κ in -3:2:3, m in -abs(κ) + 1/2:abs(κ) - 1/2
		spinor1 = SphericalSpinors(β, κ, RF(m), Pi(RF)//3, Pi(RF)//3)
		spinor2 = RadoslawSSpinors(κ, RF(m), Pi(RF)//3, Pi(RF)//3)

		spinor11 = Float64(NO(real(spinor1[1]+CF(0,0))))
		spinor12 = Float64(NO(real(spinor1[2]+CF(0,0))))

		spinor21 = Float64(NO(real(spinor2[1]+CF(0,0))))
		spinor22 = Float64(NO(real(spinor2[2]+CF(0,0))))

		@test spinor11 == spinor21
		@test spinor12 == spinor22

	end
end
###########################################################################################################
###########################################################################################################
Degree = Divide(Pi(RF),180)
x304 = real(Rotad(1,2,2,2,-1,20Degree,45Degree)) #-0.22725973883602185			# x304
x305 = real(Rotad(1,2,1,2,-1,20Degree,45Degree)) #0.14809906636301193			# x305
x306 = real(Rotad(1,2,-2,2,-1,20Degree,45Degree)) #-0.17409106008000474			# x306
x307 = real(Rotad(1,2,-2,2,-2,20Degree,45Degree)) #0.10329397779163371			# x307
x308 = real(Rotad(2,2,-2,2,-2,20Degree,45Degree)) #0.8864431717217084			# x308

@time @testset "Rotated angular functions" begin
	@test Float64(NO(real(Rotad(1,2,2,2,-1,20Degree,45Degree)))) == Float64(NO(x304))
	@test Float64(NO(real(Rotad(1,2,1,2,-1,20Degree,45Degree)))) == Float64(NO(x305))
	@test Float64(NO(real(Rotad(1,2,-2,2,-1,20Degree,45Degree)))) == Float64(NO(x306))
	@test Float64(NO(real(Rotad(1,2,-2,2,-2,20Degree,45Degree)))) == Float64(NO(x307))
	@test Float64(NO(real(Rotad(2,2,-2,2,-2,20Degree,45Degree)))) == Float64(NO(x308))
end
###########################################################################################################
###########################################################################################################
# N = n1 + n2 - 1
# Abs(l1 - l2) <= L <= l1 + l2
# -L <= M <= L
x309 = ChargeDOneC(RF(1),0,0,RF(5/10),RF(1),0,0,RF(5/10),RF(1),0,0) # 0.25						# x309
x310 = ChargeDOneC(RF(2),0,0,RF(5/10),RF(1),0,0,RF(5/10),RF(2),0,0) # 0.125						# x310
x311 = ChargeDOneC(RF(2),1,0,RF(5/10),RF(1),0,0,RF(5/10),RF(2),1,0) # 0.125						# x311
x312 = ChargeDOneC(RF(2),1,-1,RF(5/10),RF(1),1,0,RF(5/10),RF(2),2,-1) # 0.09682458365518543		# x312
x313 = ChargeDOneC(RF(2),1,-1,RF(5/10),RF(3),2,-1,RF(5/10),RF(4),1,0) # 0.05229125165837972		# x313
@time @testset "Charge density distribution" begin

@test Float64(NO(ChargeDOneC(RF(1),0,0,RF(5/10),RF(1),0,0,RF(5/10),RF(1),0,0))) == Float64(NO(x309))
@test Float64(NO(ChargeDOneC(RF(2),0,0,RF(5/10),RF(1),0,0,RF(5/10),RF(2),0,0))) == Float64(NO(x310))
@test Float64(NO(ChargeDOneC(RF(2),1,0,RF(5/10),RF(1),0,0,RF(5/10),RF(2),1,0))) == Float64(NO(x311))
@test Float64(NO(ChargeDOneC(RF(2),1,-1,RF(5/10),RF(1),1,0,RF(5/10),RF(2),2,-1))) == Float64(NO(x312))
@test Float64(NO(ChargeDOneC(RF(2),1,-1,RF(5/10),RF(3),2,-1,RF(5/10),RF(4),1,0))) == Float64(NO(x313))
end
###########################################################################################################
###########################################################################################################
x314 = OneCenterP(2,RF(2),RF(2),RF(5//10),RF(5//10),RF(20//10)) # 0.1693388246834792			# 314
x315 = OneCenterP(30,RF(20),RF(20),RF(5//10),RF(5//10),RF(20//10)) # 9.550586421859595e-34		# 315
x316 = OneCenterP(30,RF(200),RF(200),RF(5//10),RF(5//10),RF(20//10)) # 1.5356143951896117e-71	# 316
x317 = OneCenterP(100,RF(200),RF(200),RF(5//10),RF(5//10),RF(20//10)) # 4.0392168871657834e-227	# 317
x318 = OneCenterP(120,RF(200),RF(200),RF(5//10),RF(5//10),RF(20//10)) # 2.486848412057007e-270	# 318

@time @testset "One center potential" begin
	@test Float64(NO(OneCenterP(2,RF(2),RF(2),RF(5//10),RF(5//10),RF(20//10)))) == Float64(NO(x314))
	@test Float64(NO(OneCenterP(30,RF(20),RF(20),RF(5//10),RF(5//10),RF(20//10)))) == Float64(NO(x315))
	@test Float64(NO(OneCenterP(30,RF(200),RF(200),RF(5//10),RF(5//10),RF(20//10)))) == Float64(NO(x316))
	@test Float64(NO(OneCenterP(100,RF(200),RF(200),RF(5//10),RF(5//10),RF(20//10)))) == Float64(NO(x317))
	@test Float64(NO(OneCenterP(120,RF(200),RF(200),RF(5//10),RF(5//10),RF(20//10)))) == Float64(NO(x318))
end
###########################################################################################################
###########################################################################################################
x319 = Hypergeometric0F1(RF(1),RF(15//10),sprec) # 3.1655890675997798							# 319
x320 = Hypergeometric1F1(RF(1),RF(2),RF(30//10),sprec) # 6.361845641062556						# 320
x321 = HypergeometricU(RF(3),RF(2),RF(10//10)) # 0.10547895651520889							# 321
x322 = Hypergeometric2F1(CF(2),CF(3),CF(4),CF(5),sprec) # 0.15654212933375475 + 0.15079644737231007im

@time @testset "Herpergeometric functions" begin
	@test Float64(NO(Hypergeometric0F1(RF(1),RF(15//10),sprec))) == Float64(NO(x319))
	@test Float64(NO(Hypergeometric1F1(RF(1),RF(2),RF(30//10),sprec))) == Float64(NO(x320))
	@test Float64(NO(HypergeometricU(RF(3),RF(2),RF(10//10)))) == Float64(NO(x321))
	@test ComplexF64(NO(Hypergeometric2F1(CF(2),CF(3),CF(4),CF(5),sprec))) == ComplexF64(NO(x322))
end
###########################################################################################################
###########################################################################################################
x323 = AuxiliaryGk(RF(11//10),2,RF(21//10),RF(31//10),RF(41//10),RF(51//10)) # 0.2217712191193439	   x323
x324 = AuxiliaryGl(RF(11//10),2,RF(21//10),RF(31//10),RF(41//10),RF(51//10),15) # 0.010798691551480687 x324
x325 = AuxiliaryGh(RF(11//10),2,RF(21//10),RF(31//10),RF(41//10),RF(51//10),15) # 0.32834600685043897  x325
x326 = AuxiliaryG(RF(11//10),RF(21//10),RF(31//10),RF(41//10),RF(51//10),15) # 0.0277800747689916	   x326

@time @testset "Molecular sub-auxiliary functions" begin
	@test Float64(NO(
		AuxiliaryGk(RF(11//10),2,RF(21//10),RF(31//10),RF(41//10),RF(51//10)))) == Float64(NO(x323))
	@test Float64(NO(
		AuxiliaryGl(RF(11//10),2,RF(21//10),RF(31//10),RF(41//10),RF(51//10),15))) == Float64(NO(x324))
	@test Float64(NO(
		AuxiliaryGh(RF(11//10),2,RF(21//10),RF(31//10),RF(41//10),RF(51//10),15))) == Float64(NO(x325))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),RF(21//10),RF(31//10),RF(41//10),RF(51//10),15))) == Float64(NO(x326))
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
########################################### GAUX_P12_BSREP.JL #############################################
x327 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),5)# 1.3000782259762316e10  x327
x328 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),10)# -2.6937350037251043e9 x328
x329 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),15)# 6.907914778033908e6   x329
x330 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),20)# -88.57951040041706    x330
x331 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),25)# 5.1376884195114805	   x331
x332 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),30)# 5.137688371072059	   x332
x333 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),35)# 5.13768837107201	   x333
x334 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),40)#5.13768837107201	   x334
x335 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),45)#5.13768837107201	   x335
x336 = AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),50)#5.13768837107201	   x336

@time @testset "Molecular auxiliary functions gaux_p12_bsrep" begin
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),5))) == Float64(NO(x327))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),10))) == Float64(NO(x328))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),15))) == Float64(NO(x329))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),20))) == Float64(NO(x330))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),25))) == Float64(NO(x331))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),30))) == Float64(NO(x332))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),35))) == Float64(NO(x333))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),40))) == Float64(NO(x334))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),45))) == Float64(NO(x335))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),50))) == Float64(NO(x336))
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
########################################### GAUX_P123_BSREP.JL ############################################
x337 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),10)# 942.2219506111262		   x337
x338 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),15)# 968.7227043035291		   x338
x339 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),20)# 969.1689349634697		   x339
x340 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),25)# 969.169585087983		   x340
x341 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),30)# 969.169586170002		   x341
x342 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),35)# 969.1695861701844		   x342
x343 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),40)# 969.1695861701844		   x343
x344 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),45)# 969.1695861701844		   x344
x345 = AuxiliaryG(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),50)# 969.1695862			   x345

@time @testset "Molecular auxiliary functions gaux_p123_bsrep" begin
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),10))
		) == Float64(NO(x337))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),15))
		) == Float64(NO(x338))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),20))
		) == Float64(NO(x339))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),25))
		) == Float64(NO(x340))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),30))
		) == Float64(NO(x341))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),35))
		) == Float64(NO(x342))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),40))
		) == Float64(NO(x343))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),45))
		) == Float64(NO(x344))
	@test Float64(NO(
		AuxiliaryG(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),50))
		) == Float64(NO(x345))
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
############################################ GAUX_P123_REC.JL #############################################
x346 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),10)# -2.1446355854005942e36   x346
x347 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),20)# -2.793929399000712e23	   x347
x348 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),30)# 5.959231185626135e13	   x348
x349 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),40)# 845216.2070820695		   x349
x350 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),50)# 969.1706829098309		   x350
x351 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),60)# 969.1695861708716		   x351
x352 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),70)# 969.1695861701845		   x352
x353 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),80)# 969.1695861701844		   x353
x354 = AuxiliaryGr(
	RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),90)# 969.1695861701844		   x354

@time @testset "Molecular auxiliary functions gaux_p123_rec" begin
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),10))
		) == Float64(NO(x346))
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),20))
		) == Float64(NO(x347))
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),30))
		) == Float64(NO(x348))
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),40))
		) == Float64(NO(x349))
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),50))
		) == Float64(NO(x350))
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),60))
		) == Float64(NO(x351))
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),70))
		) == Float64(NO(x352))
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),80))
		) == Float64(NO(x353))
	@test Float64(NO(
		AuxiliaryGr(RF(11//10),10,RF(21//10),RF(31//10),RF(41//10),RF(51//10),RF(61//10),90))
		) == Float64(NO(x354))
end
###########################################################################################################
###########################################################################################################
###########################################################################################################
####################################### STO_MOL_INTEG_ONE_ELECT.JL ########################################
x355 = TwoCenterOverlap(
	RF(63//10),5,4,
	RF(11//2),4,4,RF(150//10),
	RF(1//10),2*Pi(RF)//3,4*Pi(RF)//3,1,15) + CF(0,1) # 0.01850589548808461							x355
x356 = TwoCenterOverlap(
	RF(36//10),2,1,
	RF(200000000000000000000001//100000000000000000000000),1,1,
	RF(8//10),RF(3//10),2*Pi(RF)//5,Pi(RF),1,15) + CF(0,1) # 0.06496217354752934					x356
x357 = TwoCenterOverlap(
	RF(300000000000000000000001//100000000000000000000000),2,1,
	RF(300000000000000000000001//100000000000000000000000),2,1,
	RF(250//10),RF(6//10),0Degree,0Degree,1,20) + CF(0,1) # -0.0004435580316765528					x357

@time @testset "Molecular integrals over Slater-type orbitals. Overlap integrals" begin
	@test Float64(real(NO(TwoCenterOverlap(
		RF(63//10),5,4,
		RF(11//2),4,4,
		RF(150//10),RF(1//10),2*Pi(RF)//3,4*Pi(RF)//3,1,15) + CF(0,1) ))
		) == Float64(real(NO(x355)))
	@test Float64(real(NO(TwoCenterOverlap(
		RF(36//10),2,1,
		RF(200000000000000000000001//100000000000000000000000),1,1,
		RF(8//10),RF(3//10),2*Pi(RF)//5,Pi(RF),1,15) + CF(0,1) ))
		) == Float64(real(NO(x356)))
	@test Float64(real(NO(TwoCenterOverlap(
		RF(300000000000000000000001//100000000000000000000000),2,1,
		RF(300000000000000000000001//100000000000000000000000),2,1,
		RF(250//10),RF(6//10),0Degree,0Degree,1,20) + CF(0,1) ))
		) == Float64(real(NO(x357)))
end
###########################################################################################################
###########################################################################################################
x358 = TwoCenterNucAttractAAB(
	RF(30000000000000000000000000001//10000000000000000000000000000),2,2,RF(124//10),
	RF(30000000000000000000000000001//10000000000000000000000000000),2,2,RF(106//10),
	RF(61//10),0Degree,0Degree,15) # 0.16031672157866095											x358
x359 = TwoCenterNucAttractAAB(
	RF(40000000000000000000000000001//10000000000000000000000000000),3,2,RF(159//10),
	RF(50000000000000000000000000001//10000000000000000000000000000),3,3,RF(107//10),
	RF(155//10),40Degree,30Degree,15) # -4.653856766826447e-6										x359
x360 = TwoCenterNucAttractAAB(
	RF(80000000000000000000000000001//10000000000000000000000000000),7,7,RF(214//10),
	RF(70000000000000000000000000001//10000000000000000000000000000),6,6,RF(208//10),
	RF(532//10),70Degree,210Degree,15) # 5.16357446188277e-5										x360
x361 = TwoCenterNucAttractAAB(
	RF(100000000000000000000000000001//10000000000000000000000000000),9,-7,RF(125//10),
	RF(100000000000000000000000000001//10000000000000000000000000000),8,-8,RF(102//10),
	RF(1007//10),80Degree,240Degree,15) # -1.5861518962960975e-6									x361

@time @testset "Molecular integrals over Slater-type orbitals. Nuclear attraction AAB integrals" begin
	@test Float64(NO(TwoCenterNucAttractAAB(
		RF(30000000000000000000000000001//10000000000000000000000000000),2,2,RF(124//10),
		RF(30000000000000000000000000001//10000000000000000000000000000),2,2,RF(106//10),
		RF(61//10),0Degree,0Degree,15))) == Float64(NO(x358))
	@test Float64(NO(TwoCenterNucAttractAAB(
		RF(40000000000000000000000000001//10000000000000000000000000000),3,2,RF(159//10),
		RF(50000000000000000000000000001//10000000000000000000000000000),3,3,RF(107//10),
		RF(155//10),40Degree,30Degree,15))) == Float64(NO(x359))
	@test Float64(NO(TwoCenterNucAttractAAB(
		RF(80000000000000000000000000001//10000000000000000000000000000),7,7,RF(214//10),
		RF(70000000000000000000000000001//10000000000000000000000000000),6,6,RF(208//10),
		RF(532//10),70Degree,210Degree,15))) == Float64(NO(x360))
	@test Float64(NO(TwoCenterNucAttractAAB(
		RF(100000000000000000000000000001//10000000000000000000000000000),9,-7,RF(125//10),
		RF(100000000000000000000000000001//10000000000000000000000000000),8,-8,RF(102//10),
		RF(1007//10),80Degree,240Degree,15))) == Float64(NO(x361))
end
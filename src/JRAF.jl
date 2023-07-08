module JRAF

using Nemo
using Cuba
using WignerSymbols
using AssociatedLegendrePolynomials
using SphericalHarmonics
using Base.Threads
using DelimitedFiles
# export JULIA_NUM_THREADS=8
sprec = 750;
ARBF = ArbField(sprec);
ACBF = AcbField(sprec);
ComplexField(prec::Int) = AcbField(prec);
CF = ComplexField(sprec);
RealField(prec::Int) = ArbField(prec);
RF = RealField(sprec);
export sprec, RF, CF, Plm, Nlm, computeYlm, SphericalHarmonics

    include("math.jl")
    include("angular_coefficients.jl")
    include("radial_coefficients.jl")
    include("special_functions.jl")
    include("gaux_p12_bsrep.jl")
    include("gaux_p123_bsrep.jl")
    include("gaux_p123_rec.jl")
    include("cgaux_p123_num.jl")
    include("sto_mol_integ_one_elect.jl")

end

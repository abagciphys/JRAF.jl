# JRAF

[![Build Status](https://travis-ci.com/abagciphys/JRAF.jl.svg?branch=master)](https://travis-ci.com/abagciphys/JRAF.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/abagciphys/JRAF.jl?svg=true)](https://ci.appveyor.com/project/abagciphys/JRAF-jl)
[![Coverage](https://coveralls.io/repos/github/abagciphys/JRAF.jl/badge.svg?branch=master)](https://coveralls.io/github/abagciphys/JRAF.jl?branch=master)

This is a package for calculation of relativistic molecular auxiliary functions derived in a paper Phys Rev E 91:023303, 2015) by Bagci and Hoggan. These auxiliary functions provide efficient computation of integrals arising at the self-consistent field level for molecules using Slater-type bases with non-integer principal quantum numbers. It applies both in relativistic and non-relativistic electronic structure theory. Highly accurate results can be achieved for these auxiliary functions via the procedure described in papers: 
<ul>
<li> Rendiconti Lincei. Scienze Fisiche e Naturali volume 29, pages191–197(2018) A. Bagci, P. E. Hoggan </li>
<li> Rendiconti Lincei. Scienze Fisiche e Naturali volume 29, pages765–775(2018) A. Bagci, P. E. Hoggam and M. Adak </li>
<li> Rendiconti Lincei. Scienze Fisiche e Naturali volume 31, pages1089–1103(2020) A. Bagci, P. E. Hoggan, <br />
and using Nemo Package developed by  William Hart, Tommy Hofmann, Claus Fieker, Fredrik Johansson. </li>
</ul>
The JRAF package inludes a "math.jl" file. In this file the basic mathematical functions are defined with the notation used by Mathematica. It is written for those familiar with the Mathematica programming language. <br />
For accuracy, both arbitrary-precision floating-point ball arithmetic which supports real and complex numbers and the MPFR library are used.

The list of special functions such as Legendre, Laguerre polynomials, spherical harmonics, spherical spinors, rotated-angular functions, hypergeometric functions ect. defined in the package can be found in "special_functions.jl" file. The mathematica notation is used for the definition of these functions too. 

The radial part of Slater-type orbitals and Slater-type spinor orbitals can be found in "radial_coefficients.jl". 

See "angular_coefficients.jl" file for the angular momentum coefficients related with product of two spherical harmonics located on different centers, the Clebsch-Gordan and Gaunt coefficients.

See "gaux_p12_bsrep.jl", "gaux_p123_bsrep.jl", "gaux_p123_rec.jl" files finally for the details on relativistic molecular auxiliary functions.

For now, in the JRAF package we focused on computation of the non-relativsitic two-center one-electron overlap, nuclear attraction and kinetic energy integrals over the real spherical harmonics. The details are in "sto_mol_integ_one_elect.jl". file.

The JRAF package will be updated periodically. 

The documentation will be available soon. 

Install with the new package manager via ]add JRAF or

## Installation
Install with the new package manager via
```julia
using Pkg
Pkg.add(path="https://github.com/abagciphys/JRAF.jl.git")
```
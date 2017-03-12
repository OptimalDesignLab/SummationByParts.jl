# SummationByParts
 
[![Build Status](https://travis-ci.org/OptimalDesignLab/SummationByParts.jl.svg?branch=master)](https://travis-ci.org/OptimalDesignLab/SummationByParts.jl)
[![codecov.io](http://codecov.io/github/OptimalDesignLab/SummationByParts.jl/coverage.svg?branch=master)](http://codecov.io/github/OptimalDesignLab/SummationByParts.jl?branch=master)

<!--- This code is for the previous, private repo
[![Build Status](https://magnum.travis-ci.com/OptimalDesignLab/SummationByParts.jl.svg?token=EpgqD9NsMEGcsBnVGzrH&branch=master)](https://magnum.travis-ci.com/OptimalDesignLab/SummationByParts.jl)
[![Coverage Status](https://coveralls.io/repos/OptimalDesignLab/SummationByParts.jl/badge.svg)](https://coveralls.io/r/OptimalDesignLab/SummationByParts.jl)
--->

## Introduction

SummationByParts is a [Julia](http://julialang.org) package that implements summation-by-parts (SBP) operators, which can be used to construct stable, high-order discretizations of partial differential equations.  SBP operators are finite-difference operators, but they share much in common with finite-element operators. For more information about SBP methods, the following reviews are a great place to start:
* M. Svärd and Jan Nordström, <a href='http://dx.doi.org/10.1016/j.jcp.2014.02.031'>"Review of summation-by-parts schemes for initial–boundary-value problems,"</a> <em>Journal of Computational Physics</em>, July, 2014.<br><br>
* D. Del Rey Fernández, J. Hicken, and D. Zingg, <a href='http://dx.doi.org/10.1016/j.compfluid.2014.02.016'>"Review of summation-by-parts operators with simultaneous approximation terms for the numerical solution of partial differential equations,"</a> <em>Computers & Fluids</em>, May, 2014.<br><br>

This package focuses on multidimensional SBP operators for the triangle and tetrahedral.  For the theory behind multidimensional SBP operators please see
* J. Hicken, D. Del Rey Fernández, and D. Zingg, <a href='http://dx.doi.org/10.1137/15m1038360'>"Multidimensional Summation-by-Parts Operators: General Theory and Application to Simplex Elements,"</a> <em>SIAM Journal on Scientific Computing</em>, July 06, 2016.<br><br>


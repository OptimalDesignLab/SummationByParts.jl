# SummationByParts

## Introduction

SummationByParts is a [Julia](http://julialang.org) package that implements summation-by-parts (SBP) operators, which can be used to construct stable, high-order discretizations of partial differential equations. SBP operators are finite-difference operators, but they share much in common with finite-element operators. For more information about SBP methods, the following reviews are a great place to start:

* M. Svärd and Jan Nordström, ["Review of summation-by-parts schemes for initial–boundary-value problems,"]('http://dx.doi.org/10.1016/j.jcp.2014.02.031') *Journal of Computational Physics*, July, 2014.

* D. Del Rey Fernández, J. Hicken, and D. Zingg, ["Review of summation-by-parts operators with simultaneous approximation terms for the numerical solution of partial differential equations,"]('http://dx.doi.org/10.1016/j.compfluid.2014.02.016') *Computers & Fluids*, May, 2014.

This package focuses on multidimensional SBP operators for the triangle and tetrahedral.  For the theory behind multidimensional SBP operators, please see:

* J. Hicken, D. Del Rey Fernández, and D. Zingg, ["Multidimensional Summation-by-Parts Operators: General Theory and Application to Simplex Elements,"]('http://dx.doi.org/10.1137/15m1038360') *SIAM Journal on Scientific Computing*, July 06, 2016.

* D. Del Rey Fernández, J. Hicken, and D. Zingg, ["Simultaneous Approximation Terms for Multi-dimensional Summation-by-Parts Operators,"]('https://doi.org/10.1007/s10915-017-0523-7') *Journal of Scientific Computing*, 2018.

SummationByParts also provides functionality to construct high-order symmetric quadrature rules with positive weights on simplices. For details on construction of such quadrature rules, please see: 

* Z. Worku, J. Hicken, D. Zingg, ["Quadrature Rules on Triangles and Tetrahedra for Multidimensional Summation-By-Parts Operators,"](https://arxiv.org/abs/2311.15576) *Submitted to Journal of Scientific Computing*, 2024.

## Using the Package

The following documentation provides a brief overview of how to use the SummationByParts package.  This assumes the user has some familiarity with Julia. 

### Building SBP operators

The construction of an SBP operator is best explained with an example.  The following code produces a degree 3 (order 4) SBP operator on a triangle.

    using SummationByParts
    sbp = SummationByParts.getTriSBPDiagE(degree=3)
    
Here is another example, which shows how to construct a degree 2 SBP element on a tetrahedron (in this example, it has been assumed that the `using SummationByParts` statement has already been executed).

    sbp = SummationByParts.getTetSBPDiagE(degree=2,Tsbp=Float64)

The SBP operators are parametrized; hence, one can specify the their type as examplified by the use of `Float64`.  This means that the matrix fields inside the `sbp` type are `Float64` arrarys. 

In general, we recommend using the SBP operators with methods provided by the package and that users *do not* rely on the fields of the sbp type directly.  This is because we may, in the future, change the fields to support different SBP operators (or some fields may become obsolete).  Nevertheless, for those who are curious, here is a partial list of the most important fields.
* `sbp.degree` : maximum polynomial degree for which the derivatives are exact
* `sbp.numnodes` : number of nodes for the operator
* `sbp.vtx` : vertices of the reference element in computational space
* `sbp.w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate 

For more information on the usage of the SummationByPart.jl package, please refer to the [examples](https://github.com/OptimalDesignLab/SummationByParts.jl/tree/master/examples).
# Reference

## Contents

```@contents
Pages = ["reference.md"]
```

## Index 

```@index
Pages = ["reference.md"]
```

## Build Operators
Functions used to build the SBP operators.
```@meta
# buildoperators.jl
```
```@docs
SummationByParts.bndrynodalexpansion
SummationByParts.nodalexpansion
SummationByParts.boundaryoperators
SummationByParts.boundarymassmatrix
SummationByParts.accuracyconstraints
SummationByParts.commuteerror
SummationByParts.buildoperators
SummationByParts.buildsparseoperators
SummationByParts.buildMinConditionOperators
SummationByParts.getnodepermutation
SummationByParts.buildoperators_pocs
```

```@meta
# buildfaceoperators.jl
```
```@docs
SummationByParts.buildfacereconstruction
SummationByParts.buildfacederivatives
```

## Cubature Module
Funtions for constructing cubature rules.
```@meta
# cubature.jl
```
```@autodocs
Modules = [SummationByParts.Cubature]
```

```@meta
# derivecubature.jl
```
```@docs
SummationByParts.deriveTriCubatureOmega
SummationByParts.deriveTriCubatureGamma
SummationByParts.deriveTriCubatureDiagE
SummationByParts.deriveTetCubatureOmega
SummationByParts.deriveTetCubatureGamma
SummationByParts.deriveTetCubatureDiagE
```

## Differentiation
Functions related to strong and weak differentiation using the SBP operators

```@meta
# differentiate
This file gathers together functions related to strong differentiation using
the SBP operators
```
```@docs
SummationByParts.differentiate!
SummationByParts.differentiateElement!
```

```@meta
# differentiate_rev
# This file contains the reverse-mode version of the methods in differentiate.jl
```
```@docs
SummationByParts.differentiate_rev!
SummationByParts.differentiateElement_rev!
```

```@meta
# directionaldifferentiate.jl
```
```@docs
SummationByParts.directionalDifferentiateElement!
```

```@meta
# weakdifferentiate.jl
# This file gathers together functions related to "weak" differentiation using
# the SBP operators
```
```@docs
SummationByParts.weakdifferentiate!
SummationByParts.weakDifferentiateElement!
```

```@meta
# weakdifferentiate_rev.jl
# This file contains the reverse-mode version of the methods in
# weakdifferentiate.jl
```
```@docs
SummationByParts.weakdifferentiate_rev!
SummationByParts.weakDifferentiateElement_rev!
```

```@meta
# weakdifferentiate_jac.jl
# This file gathers together functions related to forming Jacobians of "weak"
# differentiation using the SBP operators
```
```@docs
SummationByParts.weakDifferentiateElement_jac!
```

```@meta 
# edgestabilize.jl
# Applies edge stabilization to a given field, differentiating in the
# direction specified by `dirvec`, and scaling by the `tau` field.
```
```@docs 
SummationByParts.edgestabilize!
```

## Integration 
Functions related to volumne and face integration over a test function using SBP operators.

```@meta
# faceintegrate.jl
# This file gathers together methods related to integration over faces, both
# against test functions and for integral functionals.
```
```@docs
SummationByParts.integratefunctional!
SummationByParts.integrateBoundaryFunctional!
SummationByParts.boundaryintegrate!
SummationByParts.boundaryFaceIntegrate!
SummationByParts.interiorfaceintegrate!
SummationByParts.interiorFaceIntegrate!
```

```@meta
# faceintegrate_rev.jl
# This file contains the reverse-mode version of the methods in faceintegrate.jl
```
```@docs
SummationByParts.integratefunctional_rev!
SummationByParts.integrateBoundaryFunctional_rev!
SummationByParts.boundaryintegrate_rev!
SummationByParts.boundaryFaceIntegrate_rev!
SummationByParts.interiorfaceintegrate_rev!
SummationByParts.interiorFaceIntegrate_rev!
```

```@meta
# faceintegrate_jac.jl
# This file gathers together methods related to computing the Jacobian of
# face-based integral terms; it combines operations from both faceinterpolate
# and faceintegrate.
```
```@docs
SummationByParts.boundaryFaceIntegrate_jac!
SummationByParts.interiorFaceIntegrate_jac!
```

```@meta
# volumeintegrate.jl
# This file gathers together functions related to volumne integration over a
# test function using SBP operators
```
```@docs
SummationByParts.volumeintegrate!
SummationByParts.volumeIntegrateElement!
```

```@meta
# volumeintegrate_rev.jl
# This file contains the reverse-mode version of the methods in
# volumeintegrate.jl
```
```@docs
SummationByParts.volumeintegrate_rev!
SummationByParts.volumeIntegrateElement_rev!
```

## Mapping 
Functions related to the calculation of the coordinate mapping Jacobian on elements and faces.
```@meta
# mappingjacobian.jl
# This file gathers together functions related to the calculation of the
# coordinate mapping Jacobian on elements and faces
```
```@docs
SummationByParts.calcMappingJacobian!
SummationByParts.calcMappingJacobianElement!
SummationByParts.mappingjacobian!
```

```@meta
# mappingjacobian_rev.jl
# This file gathers together functions related to the reverse-mode differentiation of 
# mapping Jacobian
```
```@docs
SummationByParts.calcMappingJacobian_rev!
SummationByParts.mappingjacobian_rev!
```

```@meta
# facenormal.jl
# This file gathers together functions related to the calculation of scaled face
# normals for general (curvilinear) faces
```
```@docs
SummationByParts.calcFaceNormals!
SummationByParts.facenormal!
```

```@meta
# facenormal_rev.jl
# This file contains the reverse-mode version of the methods in
# facenormal.jl
```
```@docs
SummationByParts.calcFaceNormals_rev!
SummationByParts.facenormal_rev!
```

## Optimizer Module
Functions for optimization. 
```@meta
# optimizer.jl
```
```@autodocs
Modules = [SummationByParts.Optimizer]
```

## OrthoPoly Module
Functions for working with orthogonal polynomials.
```@meta
# orthopoly.jl
```
```@autodocs
Modules = [SummationByParts.OrthoPoly]
```

## Outer Constructors 
Constructors for the SBP operator classes.
```@meta 
# outerconstructor.jl
# This file gathers together outer constructors for the SBP operators
```
```@docs
SummationByParts.getLineSegSBPLobbato
SummationByParts.getLineSegSBPLegendre
SummationByParts.getTriSBPGamma
SummationByParts.getTriSBPOmega
SummationByParts.getTriSBPDiagE
SummationByParts.getTetSBPGamma
SummationByParts.getTetSBPOmega
SummationByParts.getTetSBPDiagE
SummationByParts.getLineSegFace
SummationByParts.TriFace
SummationByParts.getTriFaceForDiagE
SummationByParts.TetFace
SummationByParts.getTetFaceForDiagE
```

## SymCubatures Module
Types and methods for mapping between symmetry groups and nodes for cubatures.
on various domains
```@meta
# symcubatures.jl
```
```@autodocs
Modules = [SummationByParts.SymCubatures]
```

## Types
SBP abstract and concrete type definitions.
```@meta
# sbp_types.jl
# SBP abstract and concrete type definitions
```
```@docs
SummationByParts.AbstractSBP
SummationByParts.LineSegSBP
SummationByParts.TriSBP
SummationByParts.SparseTriSBP
SummationByParts.TetSBP
SummationByParts.SparseTetSBP
SummationByParts.AbstractFace
SummationByParts.DenseFace
SummationByParts.LineSegFace
SummationByParts.TriFace
SummationByParts.TetFace
SummationByParts.SparseFace
SummationByParts.TriSparseFace
SummationByParts.TetSparseFace
```

```@meta
# face_types.jl
# These are from ODLCommonTools; we can return to using that package when
# it is brought up-to-date
```
```@docs
SummationByParts.Boundary
SummationByParts.Interface
```

## Utilities
Functions that are not easily categorized.
```@meta
# utils.jl
# This file gathers together a hodge-podge of functions that are not easily
# categorized
```
```@docs
SummationByParts.getNumFaceNodes
SummationByParts.getnbrnodeindex
SummationByParts.calcnodes
SummationByParts.calcminnodedistance
SummationByParts.buildinterpolation
SummationByParts.permuteinterface!
SummationByParts.permuteface!
SummationByParts.basispursuit!
SummationByParts.calcSparseSolution!
SummationByParts.absMatrix!
SummationByParts.calcMatrixEigs!
SummationByParts.calcMatrixEigs_rev!
SummationByParts.conditionObj
SummationByParts.conditionObjGrad!
SummationByParts.eigenvalueObj
SummationByParts.eigenvalueObjGrad!
SummationByParts.truncErr
SummationByParts.computeConditionNumber
SummationByParts.pocs_sparse_s
SummationByParts.quadTruncErr
checkInteriorNodeLocaton
```

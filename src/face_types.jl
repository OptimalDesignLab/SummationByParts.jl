# These are from ODLCommonTools; we can return to using that package when it is
# brought up-to-date

"""
### ODLCommonTools.Boundary

Used to identify boundary faces in a finite-element grid.

**Fields**

* `element` : index of the element to which the boundary face belongs
* `face` : the face index of the boundary (local index to the element)

**Example**

To mark face 2 of element 7 to be a boundary face, use `Boundary(7,2)`

"""
immutable Boundary
  element::UInt32
  face::UInt8
end

"""
### ODLCommonTools.Interface

Used to identify interfaces between elements in a finite-element grid.

**Fields**

* `elementL` : index of the so-called left element in the pair
* `elementR` : index of the so-called right element in the pair
* `faceL` : the face index of the interface with respect to the left element
* `faceR` : the face index of the interface with respect to the right element
* `orient` : orientation of the 'right' element relative to the 'left'

**Example**

Consider an interface between elements 2 and 5.  Suppose the interface is on
face 1 of element 2 and face 3 of element 5.  Furthermore, suppose element 5 has
orientation 1 relative to element 1 (defintion of orientation TBD).  This can be
indicated as `Interface(2,5,1,3,1)`

"""
immutable Interface
  elementL::UInt32
  elementR::UInt32
  faceL::UInt8
  faceR::UInt8
  orient::UInt8
end

function show(io::IO, object::Boundary)
  print(io, "Boundary element, face = ", object.element, ", ", object.face)
end

function show(io::IO, obj::Interface)
  print(io, "Interface elementL, elementR, faceL, faceR, orient = ",
        obj.elementL, ", ",obj.elementR, ", ", obj.faceL, ", ", obj.faceR,
        ", ",obj.orient)
end

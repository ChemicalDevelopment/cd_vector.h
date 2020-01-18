# cd\_vector.h

`cd_vector.h` (aka ChemicalDevelopment Vector) is a 2d, 3d, and 4d vector/matrix library in C.


It is a single file, header only library, defining `static inline` C functions (so that in many inner loops, there is not the penalty of calling a function, saving performance in some applications).

## Types

Types beginning with `v` are vectors, i.e. `v2` is a vector of 2 floats, etc. Matrices are `vNxM`, where N is the number of rows, and M is the number of columns.

Note that these are all fixed size, so no pointer derefs are needed (nor allocation/deallocation)


  * `v2`: 2D vector of floats. All pointwise operations (`v2_add`, `v2_sub`, `v2_mul`, `v2_div`, ...) are supported
  * `v3`: 3D vector of floats. All operations that `v2` has, plus a cross product, and some 3D-specific geometric functions
  *








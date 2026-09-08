```@meta
CurrentModule = MobiusTransformations
```

# Linear algebra

A transformation carries a `2×2` coefficient matrix `[a b; c d]`, with the
usual determinant and projective normalization.

```julia
m = Möbius(1, 2, 3, 4)
Matrix(m)          # [1 2; 3 4]
det(m)             # 1*4 - 2*3
normalize(m)       # det == 1
```

```@docs
Base.Matrix
det
normalize
Base.eltype
Base.broadcastable
```

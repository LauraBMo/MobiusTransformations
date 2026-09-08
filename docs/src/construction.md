```@meta
CurrentModule = MobiusTransformations
```

# Construction

A Möbius transformation is built from its four coefficients, from a `2×2`
matrix, or from the images of three points. `Mobius` is an ASCII alias for
`Möbius`.

```julia
m = Möbius(1, 2, 3, 4)          # z -> (z + 2) / (3z + 4)
n = Möbius([1 2; 3 4])          # same map, from its matrix
f = Möbius(0, 1, Inf, 1, 2, 3)  # 0 -> 1, 1 -> 2, Inf -> 3
```

```@docs
MöbiusTransformation
Möbius
Mobius
```

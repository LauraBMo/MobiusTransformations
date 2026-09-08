```@meta
CurrentModule = MobiusTransformations
```

# Operations

Transformations compose, invert, and compare as projective maps — equal up to a
scalar multiple of the coefficients.

```julia
g = Möbius(0, 1, 1, 0)            # z -> 1/z
isone(g * inv(g))                  # true
Möbius(1, 2, 3, 4) == Möbius(2, 4, 6, 8)   # scalar multiples are equal
```

```@docs
Base.inv
Base.:(*)
Base.:∘
Base.isone
Base.:(==)
Base.hash
```

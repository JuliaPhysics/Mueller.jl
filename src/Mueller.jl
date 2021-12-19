module Mueller

using StaticArrays

export rotate, linear_polarizer, hwp, qwp, mirror, waveplate

# Various formulae for polarization components
"""
    Mueller.rotation([T=Float64], θ)

Generate a rotation matrix with the given angle, in radians, `θ`. This can be used to rotate the axes of polarization components arbitrarily. For convenience, the [`rotate`](@ref) method will rotate a component, without having to generate this matrix, itself.

# Examples

Rotate a linear polarizer by 90 degrees counter-clockwise
```jldoctest
julia> M = linear_polarizer();

julia> r = Mueller.rotation(π/4);

julia> Mr = r' * M * r;

julia> Mr ≈ linear_polarizer(π/4)
true
```

# See also
[`rotate`](@ref)
"""
function rotation(T, θ)
    sin2t, cos2t = sincos(2 * θ)
    return SA{T}[1 0 0 0
              0 cos2t sin2t 0
              0 -sin2t cos2t 0
              0 0 0 1]
end
rotation(θ) = rotation(Float64, θ)

"""
    rotate(M, θ)

Rotates the component represented by the Mueller matrix `M` counter-clockwise by angle `θ` (in radians).

# Examples
```jldoctest
julia> M = linear_polarizer();

julia> Mr = rotate(M, π/2);

julia> Mr ≈ linear_polarizer(π/2)
true
```

# See also
[`Mueller.rotation`](@ref)
"""
function rotate(mat::AbstractMatrix{T}, θ) where T
    r = rotation(T, θ)
    return r' * mat * r
end

"""
    linear_polarizer([T=Float64], θ=0; p=1)

A linear polarizer with the throughput axis given by `θ`, in radians, by default horizontal. The partial polarization can be given with the `p` keyword argument, which changes the intensity by a factor of `p^2/2`.

# Examples

```jldoctest
4×4 StaticArrays.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
  0.5   0.5   0.0  0.0
  0.5   0.5   0.0  0.0
 -0.0  -0.0  -0.0  0.0
  0.0   0.0   0.0  0.0

julia> S = [1, 0, 0, 0]; # I, Q, U, V

julia> M * S # only horizontal component (+Q) remains
4-element StaticArrays.SVector{4, Float64} with indices SOneTo(4):
 0.5
 0.5
 0.0
 0.0
```
"""
function linear_polarizer(T::Type, θ=0; p=1)
    I = T(p^2 / 2)
    sin2t, cos2t = sincos(2 * θ)
    M = I * SA{T}[1 cos2t sin2t 0
                 cos2t cos2t^2 cos2t * sin2t 0
                 -sin2t -cos2t * sin2t -sin2t^2 0
                 0 0 0 0]
    return M
end
linear_polarizer(θ=0; kwargs...) = linear_polarizer(Float64, θ; kwargs...)

# wave plates

"""
    waveplate([T=Float64], θ=0, δ=0; p=1)

A generic phase retarder (waveplate) with fast axis aligned with angle `θ`, in radians, and phase delay of `δ`, in radians, along the slow axis. The degree of polarization can be set with the keyword argument `p`, which will change the intensity by a factor of `p^2/2`.

# Examples
```jldoctest
julia> M = waveplate(0, π);

julia> M ≈ hwp()
true

julia> M = waveplate(π/2, π/2);

julia> M ≈ qwp(π/2)
true
```

# See also
[`hwp`](@ref), [`qwp`](@ref), [`mirror`](@ref)
"""
function waveplate(T::Type, θ=0, δ=0; p=1)
    sin2t, cos2t = sincos(2 * θ)
    sind, cosd = sincos(δ)
    I = T(p^2 / 2)
    return I * SA{T}[1 0 0 0
                     0 cos2t^2 + sin2t^2*cosd cos2t * sin2t * (1 - cosd) sin2t * sind
                     0 cos2t*sin2t*(1-cosd) cos2t^2*cosd + sin2t^2 -cos2t * sind
                     0 -sin2t * sind cos2t * sind cosd]
end
waveplate(θ=0, δ=0; kwargs...) = waveplate(Float64, θ, δ; kwargs...)

"""
    hwp([T=Float64], θ=0; p=1)

A half-wave plate (HWP) with fast axis oriented at angle `θ`, in radians. The degree of polarization can be set with the keyword argument `p`, which will change the intensity by a factor of `p^2/2`.

# Examples
```jldoctest
julia> hwp()
4×4 StaticArrays.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 0.5   0.0   0.0           0.0
 0.0   0.5   0.0           0.0
 0.0   0.0  -0.5          -6.12323e-17
 0.0  -0.0   6.12323e-17  -0.5

julia> S = [1, 1, 0, 0]; # I, Q, U, V

julia> M * S # allow +Q through unchanged
4-element StaticArrays.SVector{4, Float64} with indices SOneTo(4):
 0.5
 0.5
 0.0
 0.0

julia> rotate(M, π/8) * S # switch +Q to +U
4-element StaticArrays.SVector{4, Float64} with indices SOneTo(4):
  0.5
  9.983673087713696e-17
  0.5
 -4.329780281177466e-17

```

# See also
[`waveplate`](@ref), [`qwp`](@ref)
"""
hwp(T::Type, θ=0; kwargs...) = waveplate(T, θ, π; kwargs...)
hwp(θ=0; kwargs...) = hwp(Float64, θ; kwargs...)


"""
    qwp([T=Float64], θ=0; p=1)

A quarter-wave plate (QWP) with fast axis oriented at angle `θ`, in radians. The degree of polarization can be set with the keyword argument `p`, which will change the intensity by a factor of `p^2/2`.

# Examples
```jldoctest
4×4 StaticArrays.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 0.5   0.0   0.0           0.0
 0.0   0.5   0.0           0.0
 0.0   0.0  -0.5          -6.12323e-17
 0.0  -0.0   6.12323e-17  -0.5

julia> S = [1, 1, 0, 0]; # I, Q, U, V

julia> M * S # allow +Q through unchanged
4-element StaticArrays.SVector{4, Float64} with indices SOneTo(4):
 0.5
 0.5
 0.0
 0.0

julia> hwp(π/8) * S # switch +Q to +U
4-element StaticArrays.SVector{4, Float64} with indices SOneTo(4):
  0.5
  1.1102230246251565e-16
  0.5
 -4.329780281177466e-17
```

# See also
[`waveplate`](@ref), [`hwp`](@ref)
"""
qwp(T::Type, θ=0; kwargs...) = waveplate(T, θ, π/2; kwargs...)
qwp(θ=0; kwargs...) = qwp(Float64, θ; kwargs...)

"""
    mirror([T=Float64], r=1, θ=0, δ=π)

A reflective mirror with reflectance `r`, oriented at angle `θ`, in radians, compared to the reference frame of the light, and with phase shift `δ`. An ideal mirror will have perfect reflectance and a 180° phase shift.

# Examples
```jldoctest
julia> M = mirror()

julia> S = [1, 1, 0, 0]; # I, Q, U, V

julia> M * S # no change
4-element StaticArrays.SVector{4, Float64} with indices SOneTo(4):
 1.0
 1.0
 0.0
 0.0

julia> mirror(1, π/4) * S # rotates polarized light
4-element StaticArrays.SVector{4, Float64} with indices SOneTo(4):
  1.0
 -1.0
  1.2246467991473532e-16
  1.2246467991473532e-16
```
"""
function mirror(T::Type, r=1, θ=0, δ=π)
    a = 0.5 * (r + 1)
    b = 0.5 * (r - 1)
    sin2t, cos2t = sincos(2 * θ)
    sin4t = sin(4 * θ)
    sind, cosd = sincos(δ)
    sqrm = sqrt(r)
    M = SA{T}[a          b * cos2t                            b * sin2t                           0
              b * cos2t  a * cos2t^2 + sqrm * cosd * sin2t^2 (a - sqrm * cosd) * sin4t * 0.5     -sqrm * sind * sin2t
              b * sin2t (a - sqrm * cosd) * sin4t * 0.5       a * sin2t^2 + sqrm * cosd * cos2t^2 sqrm * sind * cos2t
              0          sqrm * sind * sin2t                 -sqrm * sind * cos2t                 sqrm * cosd]
    return M
end
mirror(r=1, θ=0, δ=π) = mirror(Float64, r, θ, δ)

end # module

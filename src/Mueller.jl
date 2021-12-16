module Mueller

using StaticArrays

export rotation, linear_polarizer

# Various formulae for polarization components
function rotation(T, θ)
    sin2t, cos2t = sincos(2 * θ)
    return SA{T}[1 0 0 0
              0 cos2t sin2t 0
              0 -sin2t cos2t 0
              0 0 0 1]
end
rotation(θ) = rotation(Float64, θ)

function linear_polarizer(T::Type, gamma=0, p=1)
    sin2g, cos2g = sincos(2 * gamma)
    I = T(p^2 / 2)
    return I * SA{T}[1 cos2g 0 0
                 cos2g 1 0 0
                 0 0 sin2g 0 
                 0 0 0 sin2g]    
end
linear_polarizer(args...) = linear_polarizer(Float64, args...)

end # module

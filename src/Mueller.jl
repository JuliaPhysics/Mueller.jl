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

function rotate(mat::AbstractMatrix{T}, θ) where T
    r = rotation(T, θ)
    return r * mat * r
end

function linear_polarizer(T::Type, gamma=0, p=1)
    I = T(p^2 / 2)
    M = I * SA{T}[1 1 0 0
                 1 1 0 0
                 0 0 0 0 
                 0 0 0 0]    
    return rotate(M, gamma)
end
linear_polarizer(gamma=0, p=1) = linear_polarizer(Float64, gamma, p)

end # module

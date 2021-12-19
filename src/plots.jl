using LinearAlgebra
using RecipesBase

# polarization ellipse
@userplot PolEllipse

"""
    polellipse(S)
    polellipse!(S)

Plot the polarizatio ellipse for the given Stokes parameters. The polarization ellipse draws the path of the eletric field over one wavelength. The input vector should be a Stokes vector `[I, Q, U, V]`. This is a user recipe and can be composed with other plot recipes.
"""
polellipse

@recipe function f(pe::PolEllipse, N=1000)
    aspect_ratio --> 1
    xguide --> "x"
    yguide --> "y"

    # extract elliptical params (a, b, θ) from Stokes vector
    stokes = first(pe.args)
    L = hypot(stokes[begin + 1], stokes[begin + 2])
    Ip = hypot(L, stokes[begin + 3])
    θ = 0.5 * atan(stokes[begin + 2], stokes[begin + 1])
    a = sqrt(0.5 * (Ip + L))
    b = sqrt(0.5 * (Ip - L))
    # draw ellipse
    t = range(0, 2π, N)
    x = @. a * cos(t) * cos(θ) - b * sin(t) * sin(θ)
    y = @. a * cos(t) * sin(θ) + b * sin(t) * cos(θ)
    xlim --> (-L, L)
    ylim --> (-L, L)
    return x, y
end
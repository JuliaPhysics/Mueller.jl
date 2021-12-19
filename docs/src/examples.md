# Examples

The following are some examples using Mueller.jl. These examples are simple demonstrations of the wave-nature of light, and how it can be measured and changed.

## Bell's inequality

Inspired by [this video](https://www.youtube.com/watch?v=zcqZHYo7ONs), we can demonstrate [Bell's Theroem](https://en.wikipedia.org/wiki/Bell%27s_theorem) on orthogonal states of polarized light. First, let's look at the effect a single linear polarizer has on unpolarized light.

```@example bell
using Mueller
using Plots
using Unitful: °
using UnitfulRecipes

# Stokes vector: I, Q, U, V
S = [1, 0, 0, 0]
M0 = linear_polarizer()
Sp = M0 * S
```

the total intensity is halved, and that half is polarized in the +Q direction. We can visualize this by drawing the polarization ellipse, which is the path traced by the electric field in the electromagnetic wave over one cycle. The [`polellipse`](@ref) plot recipe plots this ellipse from a Stokes vector.

```@example bell
polellipse(Sp, label="0° LP", lims=(-1, 1))
```

If we add another polarizer in the orthogonal direction, we will end up attenuating all of the light (by the definition of orthogonal). In order to combine Mueller matrices, we just use matrix multiplication from right to left

```math
M = M_N \cdot \left(\hdots \cdot M_2 \cdot \left(M_1 \right)\right)
```

```@example bell
M2 = linear_polarizer(90°)
M = M2 * M0
Sp = M * S
```

the intensity is 0, so there is no remaining light.

What do you expect to happen, though, if we insert a polarizer in-between the two orthogonal directions at an intermediate angle, say: 45°?

```@example bell
M1 = linear_polarizer(45°)
M = M2 * M1 * M0
Sp = M * S
```

interestingly, despite light having to pass through two orthogonal polarization states (0° and 90°), due to the probabilistic nature of light 1/8 of the intensity passes through. To get a better intuition, let's look at the polarization ellipses at each stage

```@example bell
polellipse(M0 * S, label="0° LP", lims=(-1, 1))
polellipse!(M1 * M0 * S, label="0°+45° LP")
polellipse!(Sp, label="0°+45°+90° LP")
```

TODO

```@example bell
angles = range(0°, 90°, length=100)
intens = map(angles) do θ
    M = M2 * linear_polarizer(θ) * M0
    Sp = M * S
    return Sp[1]
end
plot(angles, intens, leg=false, xlabel="θ", ylabel="I")
vline!([45°], c=:black, ls=:dash, alpha=0.7)
```

## Differential polarimetry

## Generating circularly polarized light


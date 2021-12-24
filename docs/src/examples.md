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
M = M_N \cdot \left(\ldots \cdot M_2 \cdot \left( M_1 \right)\right)
```

!!! tip "Ordering"
    Let's say you had a list of optical components in order of the light's path. To quickly produce a Mueller matrix from their combination, do

    ```julia
    julia> components = # M0, M1, M2, ...

    julia> M = prod(reverse(components))
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

We can plot the light transmitted as a function of the angle of the  intermediate polarizer

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

this curve follows the form ``\sin(2\theta)`` and is part of the proof  of  Bell's theorem (Bell's inequality) which states there are no "hidden" variables that can know the outcome of the light's polarized state before arriving at the polarizers.

## Differential polarimetry

Stokes parameters are a convenient basis for polarization due to its direct relation to observables. To demonstrate this, let's set up an experiment with a simple polarimeter. This polarimeter consists of a half-wave plate (HWP) and a linear polarizer.

```@example pdi
using Mueller
M = linear_polarizer() * hwp()
```

polarized light will enter the detector, and in the current configuration the +Q polarization will be measured

```@example pdi
S = [1, 0.3, 0.2, 0]
Sp = M * S
```

now, if we rotate the HWP by ``\gamma``, the angle of linear polarization will increase by ``2\gamma``. We can measure orthogonal polarization states (e.g., +Q and -Q) by choosing the appropriate angles for ``\gamma``. In this case, we measure -Q by rotating the HWP by 45°.

```@example pdi
using Unitful: °
M45 = linear_polarizer() * hwp(45°)
Sp45 = M45 * S
```

if we take the difference of the two measurements, `Sp` and `Sp45`, we remove the unpolarized intensity component of the light and retain only the polarized component (in this case, stokes Q)

```@example pdi
diffQ = Sp - Sp45
```

from this difference, we get a clean observable of the polarimetric signal

```@example pdi
Qhat = diffQ[2]
Qhat ≈ S[2]
```

Let's repeat this process with the HWP at 22.5° and 67.5° we can measure the +U and -U states

```@example pdi
M225 = linear_polarizer() * hwp(22.5°)
Sp225 = M225 * S
M675 = linear_polarizer() * hwp(67.5°)
Sp675 = M675 * S
diffU = Sp225 - Sp675
Uhat = diffU[2]
```

so, from four measurements, we can generate the I, Q, and U stokes parameters

```@example pdi
Ihat = 0.5 * (Sp[1] + Sp45[1] + Sp225[1] + Sp675[1])
Shat = [Ihat, Qhat, Uhat, NaN]
```

which faithfully captures the linear components of the original light

```@example pdi
Shat[1:3] ≈ S[1:3]
```

## Generating circularly polarized light

Circularly polarized light can be generated from unpolarized light with a polarizer and a quarter-wave plate (QWP).

```@example circular
using Mueller
using Unitful: °

M = qwp() * linear_polarizer(45°)
```

```@example circular
S = [1, 0, 0, 0]
Sp = M * S
```

The light is  attenuated due to the linear polarizer, which is required to apply the phase change with the QWP to appropriately create  circular polarization. We can gain more intuition for this process by looking at the polarization ellipses of the light between the components

```@example circular
using Plots

polellipse(linear_polarizer(45°) * S, label="45° LP", lims=(-1, 1))
polellipse!(Sp, label="45° LP + QWP")
```

## Symbolic calculations

Here we show the flexibility of Julia's multiple dispatch by combining [Symbolics.jl](https://github.com/SciML/Symbolics.jl) with Mueller.jl to derive an equation similar to eq. 13 from ["Mueller matrix for imperfect, rotated mirrors"](https://psfcsv10.psfc.mit.edu/~sscott/MSEmemos/mse_memo_20c.pdf)

```@example mirror
using Mueller
using Symbolics
using Unitful: °

@variables I γ χ

# represent light using the polarization ellipse angle γ
S = [I, cos(2γ), sin(2γ), 0]
# need to specify eltype for symbolic variable
M = mirror(typeof(I), 1, 180°, χ)
```

```@example mirror
Iout = M * S
```

applying the substitutions

```@example mirror
rules = [
    cos(2χ)^2 - sin(2χ)^2 => cos(4χ),
    sin(2χ)^2 - cos(2χ)^2 => sin(4χ)
]
```

shows that

```@example mirror
substitute(Iout, rules)
```

which fully recreates equation 13 from the manuscript.

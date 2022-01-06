# Examples

The following are some examples using Mueller.jl. These examples are simple demonstrations of the wave nature of light, and how it can be measured and changed.

## Bell's inequality

Inspired by [this video](https://www.youtube.com/watch?v=zcqZHYo7ONs), we can observe the effect of [Bell's Inequality](https://en.wikipedia.org/wiki/Bell%27s_theorem) on orthogonal states of polarized light. First, let's look at the effect a single linear polarizer has on unpolarized light.

```@example bell
using Mueller

# Stokes vector: I, Q, U, V
S = [1, 0, 0, 0]
M0 = linear_polarizer()
Sp = M0 * S
```

the total intensity is halved, and that half is polarized in the +Q direction. We can visualize this by drawing the polarization ellipse, which is the path traced by the electric field in the electromagnetic wave over one cycle. The [`polellipse`](@ref) plot recipe plots this ellipse from a Stokes vector.

```@example bell
using Plots
using UnitfulRecipes

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

now let's combine the two polarizors

```@example bell
using Unitful: °
# use Unitful degrees, otherwise specify in radians
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

Interestingly, despite light having to pass through two orthogonal polarization states (0° and 90°), due to the probabilistic nature of light, 1/8 of the intensity passes through. To get a better intuition for what's happening, let's look at the polarization ellipses at each step

```@example bell
polellipse(M0 * S, label="0° LP", lims=(-1, 1))
polellipse!(M1 * M0 * S, label="0°+45° LP")
polellipse!(Sp, label="0°+45°+90° LP")
```

We can rotate the intermediate polarizor and see how that changes the amount of transmitted light.

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

This curve follows the form ``\sin(2\theta)`` and is part of the proof of Bell's theorem (Bell's inequality) which states there are no "hidden" variables that can know the outcome of the light's *orthogonal* polarized state before arriving at the polarizers.

## Differential polarimetry

Stokes parameters are a convenient basis for polarization due to their direct relation to observables. To demonstrate this, let's set up an experiment with a simple polarimeter. This polarimeter consists of a half-wave plate (HWP) and a linear polarizer.

```@example pdi
using Mueller
M = linear_polarizer() * hwp()
```

Polarized light will enter the detector and, in the current configuration, the detect the +Q polarization amplitude.

```@example pdi
S = [1, 0.3, 0.2, 0]
Sp = M * S
```

Now, if we rotate the HWP by ``\gamma``, the angle of linear polarization will increase by ``2\gamma``. We can measure orthogonal polarization states (e.g., +Q and -Q) by choosing the appropriate angles for ``\gamma`` such that the angle of linear polarization increases by 90°. Therefore, we measure -Q by rotating the HWP by 45°.

```@example pdi
using Unitful: °
M45 = linear_polarizer() * hwp(45°)
Sp45 = M45 * S
```

If we take the difference of the two measurements, `Sp` and `Sp45`, we remove the unpolarized intensity component of the light and retain only the polarized component (in this case, stokes Q).

```@example pdi
diffQ = Sp - Sp45
```

From this difference, we get a clean observable of the polarimetric signal.

```@example pdi
Qhat = diffQ[2]
Qhat ≈ S[2]
```

Let's repeat this process with the HWP at 22.5° and 67.5°, which allows us to measure the +U and -U amplitudes.

```@example pdi
M225 = linear_polarizer() * hwp(22.5°)
Sp225 = M225 * S
M675 = linear_polarizer() * hwp(67.5°)
Sp675 = M675 * S
diffU = Sp225 - Sp675
Uhat = diffU[2]
```

So, from four measurements, we can generate clean I, Q, and U Stokes parameters:

```@example pdi
Ihat = 0.5 * (Sp[1] + Sp45[1] + Sp225[1] + Sp675[1])
Shat = [Ihat, Qhat, Uhat, NaN]
```

Which faithfully captures the linear components of the original light!

```@example pdi
Shat[1:3] ≈ S[1:3]
```

!!! tip "An alternative approach"
    Mathematically, we can express the entire differential approach in a single matrix, since subtraction is distributive over matrices and vectors. For example, to get `Qhat` directly-
    ```math
    \hat{Q} = M_{0} \cdot S - M_{45} \cdot S = (M_{0} - M_{45}) \cdot S
    ```

```@example pdi
Qhat_dir =  (M - M45) * S
```

You'll quickly realize, this Mueller matrix is just indexing into the Stokes vector

```@example pdi
M_Qhat = M - M45
```

which makes it somewhat trivial in terms of practical applications.

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

The light is attenuated due to the linear polarizer, which is required for the QWP to apply the phase shift which creates circular polarization. We can gain more intuition for this process by looking at the polarization ellipses of the light between each step

```@example circular
using Plots

polellipse(linear_polarizer(45°) * S, label="45° LP", lims=(-1, 1))
polellipse!(Sp, label="45° LP + QWP")
```

## Symbolic calculations

Here we show the flexibility of Julia's multiple dispatch by combining [Symbolics.jl](https://github.com/SciML/Symbolics.jl) with Mueller.jl to derive an equation similar to Eq. 13 from ["Mueller matrix for imperfect, rotated mirrors"](https://psfcsv10.psfc.mit.edu/~sscott/MSEmemos/mse_memo_20c.pdf)

```@example mirror
using Latexify
using Mueller
using Symbolics
using Unitful: °

@variables I γ χ

# represent light using the polarization ellipse angle γ
S = [I, cos(2γ), sin(2γ), 0]
# need to specify eltype for symbolic variable
# perfect mirror
M = mirror(typeof(I), 1, 180°, χ)
nothing # hide
```

```math
\begin{equation}
\left[
\begin{array}{cccc}
1.0 & 0.0 & 0.0 & 0 \\
0.0 &  - \sin^{2}\left( 2 \chi \right) + \cos^{2}\left( 2 \chi \right) & \sin\left( 4 \chi \right) & -0.0 \\
0.0 & \sin\left( 4 \chi \right) &  - \cos^{2}\left( 2 \chi \right) + \sin^{2}\left( 2 \chi \right) & 0.0 \\
0 & 0.0 & -0.0 & -1.0 \\
\end{array}
\right]
\end{equation}
```

```@example mirror
# "apply" mirror to light
Iout = M * S
nothing # hide
```

```math
\begin{equation}
\left[
\begin{array}{c}
I \\
\left(  - \sin^{2}\left( 2 \chi \right) + \cos^{2}\left( 2 \chi \right) \right) \cos\left( 2 \gamma \right) + \sin\left( 2 \gamma \right) \sin\left( 4 \chi \right) \\
\left(  - \cos^{2}\left( 2 \chi \right) + \sin^{2}\left( 2 \chi \right) \right) \sin\left( 2 \gamma \right) + \cos\left( 2 \gamma \right) \sin\left( 4 \chi \right) \\
0.0 \\
\end{array}
\right]
\end{equation}
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
Isub = substitute.(Iout, (rules,))
nothing # hide
```

```math
\begin{equation}
\left[
\begin{array}{c}
I \\
\cos\left( 2 \gamma \right) \cos\left( 4 \chi \right) + \sin\left( 2 \gamma \right) \sin\left( 4 \chi \right) \\
\cos\left( 2 \gamma \right) \sin\left( 4 \chi \right) + \sin\left( 2 \gamma \right) \sin\left( 4 \chi \right) \\
0.0 \\
\end{array}
\right]
\end{equation}
```

which fully recreates Eq. 13 from the manuscript.

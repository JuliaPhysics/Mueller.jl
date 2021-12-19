using LinearAlgebra
using Mueller
using RecipesBase: apply_recipe
using StableRNGs
using Test
using Unitful: °

rng = StableRNG(128584)

@testset "Mueller.jl" begin

    @testset "rotations" for EL in [linear_polarizer, waveplate, qwp, hwp, mirror]
        M = EL()
        for θ in 2 * π .* rand(rng, 1000)
            R = Mueller.rotation(θ)
            @test R' * M * R ≈ rotate(M, θ)
            @test rotate(rotate(M, θ), -θ) ≈ M
        end
    end

    @testset "linear polarizers" for T in [Float32, Float64, BigFloat]
        @testset "horizontal" begin
            horizontal = linear_polarizer(T)
            @test horizontal ≈ 0.5 * T[1 1 0 0;
                                       1 1 0 0
                                       0 0 0 0
                                       0 0 0 0]

            # parallal light perfect
            S = T[1, 1, 0, 0]
            Sp = horizontal * S
            @test Sp ≈ T[1, 1, 0, 0]

            # orthogonal light blocked
            S = T[1, -1, 0, 0]
            Sp = horizontal * S
            @test Sp ≈ zeros(T, 4)

            # mixed light attenuated
            S = T[1, 1/√2, 1/√2, 0]
            Sp = horizontal * S
            @test Sp ≈ 1/2 * T[1 + 1/√2, 1  + 1/√2, 0, 0]
        end
        @testset "vertical" begin
            vertical = linear_polarizer(T, 90°)
            @test vertical ≈ 0.5 * T[1 -1 0 0;
                                     -1  1 0 0
                                     0 0 0 0
                                     0 0 0 0]

            # parallal light perfect
            S = T[1, -1, 0, 0]
            Sp = vertical * S
            @test Sp ≈ T[1, -1, 0, 0]

            # orthogonal light blocked
            S = T[1, 1, 0, 0]
            Sp = vertical * S
            @test Sp ≈ zeros(T, 4)

            # mixed light attenuated
            S = T[1, -1/√2, 1/√2, 0]
            Sp = vertical * S
            @test Sp ≈ 1/2 * T[1 + 1/√2, -1 - 1/√2, 0, 0]
        end
        @testset "double polarization" begin
            M0 = linear_polarizer(T)
            M90 = linear_polarizer(T, 90°)
            M = M0 * M90
            @test all(≈(0), M)
    
            # more randomly
            for θ in 180 .* rand(rng, 1000)
                M1 = linear_polarizer(T, θ*°)
                M2 = linear_polarizer(T, (θ + 90)*°)
                M = M0 * M90
                @test all(≈(0), M)
            end
        end
    end

    @testset "waveplates" for T in [Float32, Float64, BigFloat]
        M = waveplate(T)
        @test eltype(M) == T


        @testset "hwp" begin
            M = hwp(T)
            @test eltype(M) == T
            @test M ≈ waveplate(T, 0, π)

            Sp = M * T[1, 1, 0, 0]
            @test Sp ≈ T[1, 1, 0, 0] atol=1e-10
            Sp = M * T[1, 0, 1, 0]
            @test Sp ≈ T[1, 0, -1, 0] atol=1e-10
        end

        @testset "qwp" begin
            M = qwp(T)
            @test eltype(M) == T
            @test M ≈ waveplate(T, 0, π/2)

            Sp = M * T[1, 1, 0, 0]
            @test Sp ≈ T[1, 1, 0, 0] atol=1e-10
            Sp = M * T[1, 0, 1, 0]
            @test Sp ≈ T[1, 0, 0, 1] atol=1e-10
        end

        @testset "generate circular polarization" begin
            M = qwp(T) * linear_polarizer(T, 45°)
            S = T[1, 0, 0, 0]
            Sp = M * S
            @test Sp ≈ T[0.5, 0, 0, -0.5] atol=1e-10
        end

    end

    @testset "plot recipe" begin
        S = [1, 0, 0, 1]
        obj = Mueller.PolEllipse([S])
        recipes = apply_recipe(Dict{Symbol,Any}(), obj)
        recipe = only(recipes)
        x, y = recipe.args[1:2]
        @test minimum(x) ≈ -1/sqrt(2) atol=1e-5
        @test maximum(x) ≈ 1/sqrt(2) atol=1e-5
        @test minimum(y) ≈ -1/sqrt(2) atol=1e-5
        @test maximum(y) ≈ 1/sqrt(2) atol=1e-5
    end
end

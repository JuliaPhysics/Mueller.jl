using LinearAlgebra
using Mueller
using StableRNGs
using Test
using Unitful: °

rng = StableRNG(128584)

@testset "Mueller.jl" begin

    @testset "linear polarizers" for T in [Float32, Float64, BigFloat]
        @testset "horizontal" begin
            horizontal = linear_polarizer(T)

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
    end

    @testset "double polarization" begin
        M0 = linear_polarizer()
        M90 = linear_polarizer(90°)
        M = M0 * M90
        @test all(≈(0), M)

        # more randomly
        for θ in 180 .* rand(rng, 1000)
            M1 = linear_polarizer(θ*°)
            M2 = linear_polarizer((θ + 90)*°)
            M = M0 * M90
            @test all(≈(0), M)
        end
    end
end

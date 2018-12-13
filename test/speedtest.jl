#!/usr/bin/env julia


module test_highell

using TwoFAST
using TwoFAST.Miller
using BenchmarkTools
using LinearAlgebra
using Test

using DelimitedFiles


function test_calc_f0()
    ℓ = 10_000
    R = 0.9
    Δℓ = 4
    kmin = 1e-5
    χ0 = 1e-3
    α = kmin * χ0
    q = 1.1
    G = 18.4
    m = 100
    t = 2π * m / G
    n = q - 1 - im * t
    a = n / 2 + Δℓ / 2
    BCfn(ell) = TwoFAST.BCDEfn(ℓ, Δℓ, a, R)
    fnull = [Complex{Float64}(0), Complex{Float64}(0)]

    f210 = TwoFAST.calc_f0(R, n, Δℓ)
    #writedlm("f0", f0)
    B0 = TwoFAST.Mellell_pre(0, 0+Δℓ, R, n, α)
    f0 = B0 * f210
    laminf1 = R
    laminf2 = 1 / R
    fasymp = [Complex{Float64}(1), Complex{Float64}(1)]
    fell, ell = miller(ℓ, BCfn, f0, fasymp, laminf1, laminf2)

    # test f0
    f210true = [4.3909437819996155e-5 + 9.125604981761485e-6im,
                -4.419845719133748e-6 + 1.4156528124157448e-5im]
    @show f210[1], f210[2]
    @test abs(real(f210[1] - f210true[1])) < 1e-10norm(f210true)
    @test abs(imag(f210[1] - f210true[1])) < 1e-10norm(f210true)
    @test abs(real(f210[2] - f210true[2])) < 1e-10norm(f210true)
    @test abs(imag(f210[2] - f210true[2])) < 1e-10norm(f210true)

    # test f21
    felltrue = [2.6260385192e-296 - 4.0387685923e-296im,
                2.5902318450e-296 - 4.0681095701e-296im]
    @show fell[1], fell[2], ell
    @test ell == 3171
    @test abs(real(fell[1] - felltrue[1])) < 1e-10norm(felltrue)
    @test abs(imag(fell[1] - felltrue[1])) < 1e-10norm(felltrue)
    @test abs(real(fell[2] - felltrue[2])) < 1e-10norm(felltrue)
    @test abs(imag(fell[2] - felltrue[2])) < 1e-10norm(felltrue)

    # benchmark
    @btime TwoFAST.calc_f0($R, $n, $Δℓ)
    @btime miller($ℓ, $BCfn, $f0, $fasymp, $laminf1, $laminf2)
end


function test_calc_f21_RqmG()
    ℓ = 10_000
    R = 0.9
    Δℓ = 4
    kmin = 1e-5
    χ0 = 1e-3
    α = kmin * χ0

    @btime TwoFAST.calc_2f1_RqmG($ℓ, $R, $Δℓ; q=1.1, m=100, G=18.4, alpha=$α)

    # test correctness (assuming the result was correc, and that ℓmax didn't change)
    f21, ℓmax = TwoFAST.calc_2f1_RqmG(ℓ, R, Δℓ; q=1.1, m=100, G=18.4, alpha=α)
    f21true = [2.6260385192e-296 - 4.0387685923e-296im,
               2.5902318450e-296 - 4.0681095701e-296im]
    ℓmaxtrue = 3171
    @show f21[1], f21[2], ℓmax
    @test ℓmax == ℓmaxtrue
    @test abs(real(f21[1] - f21true[1])) < 1e-10norm(f21true)
    @test abs(imag(f21[1] - f21true[1])) < 1e-10norm(f21true)
    @test abs(real(f21[2] - f21true[2])) < 1e-10norm(f21true)
    @test abs(imag(f21[2] - f21true[2])) < 1e-10norm(f21true)
end


function high_level(ℓmax=10_000)
    #RR = [0.1:0.1:0.9; [0.99, 0.999, 0.9999, 1.0, 1.1]]
    RR = [0.8]
    F21EllCache(ℓmax+2, RR, 1600; q=1.1, kmin=1e-5, kmax=1e3, χ0=1e-3, Δℓ=4, ΔℓRg1=-4)
end


end

# TODO:
# - What's the difference between test_calc_f0() and test_calc_f21_RqmG()? Make
#   a test for the change!


#test_highell.test_calc_f0()
test_highell.test_calc_f21_RqmG()
#test_highell.high_level()


# vim: set sw=4 et sts=4 :

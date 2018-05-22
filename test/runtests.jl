#!/usr/bin/env julia

include("PkSpectra.jl")
include("TwoFASTTestLib.jl")


module TestTwoFAST

using TwoFAST
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end
using Dierckx
using PkSpectra
using TwoFASTTestLib


############### test xicalc #######################

function test_xiln(ℓ)
    νν = -2:3
    r0, ξ0 = get_quadosc_xi(ℓ)
    r1, ξ1 = calc_2fast_xi(ℓ, νν)
    for i=1:length(νν)
        if !isfinite(ξ0[1,i])
            continue
        end
        ν = νν[i]
        s0 = Spline1D(r0, ξ0[:,i])
        s1 = Spline1D(r1, ξ1[i])
        r = r0[1:10:end]
        atol = 1e-5
        rtol = 1e-4
        print("Testing ξ(r), ℓ=$ℓ, ν=$ν ... ")
        @test all(isapprox.(s0(r), s1(r), atol=atol, rtol=rtol))
        println("passed")
    end
end


function test_xi_derivs()
    # get quadosc curve
    rq, ξq0 = get_quadosc_xi(0)
    ξq = Spline1D(rq, ξq0[:,3], k=5)
    ξq′(r) = derivative(ξq, r; nu=1)
    ξq″(r) = derivative(ξq, r; nu=2)

    # get 2-FAST curve
    r0, ξ00 = calc_2fast_xi(0, [0])
    r1, ξ1m1 = calc_2fast_xi(1, [-1])
    r2, ξ2m2 = calc_2fast_xi(2, [-2])
    ξ = Spline1D(r0, ξ00[1])
    ξ′ = Spline1D(r1, - ξ1m1[1] ./ r1)
    ξ″ = Spline1D(r2, (ξ2m2[1] - ξ1m1[1]) ./ r2.^2)

    # ξ
    print("Testing ξ(r) ... ")
    r = rq[1:10:end]
    @test all(isapprox.(ξ(r), ξq(r), atol=1e-5, rtol=1e-5))
    println("passed")

    # ξ′
    print("Testing ξ′(r) ... ")
    r = r[r .>= 5]
    @test all(isapprox.(ξ′(r), ξq′(r), atol=1e-5, rtol=1e-5))
    println("passed")

    # ξ″
    print("Testing ξ″(r) ... ")
    @test all(isapprox.(ξ″(r), ξq″(r), atol=1e-5, rtol=1e-5))
    println("passed")
end


function test_xi()
    for ℓ=0:4
        test_xiln(ℓ)
    end
    test_xi_derivs()
end


##################### test wljj ##############################

function test_wl_χ2303_R(R)
    d = readdlm("data/wljj_chi2303.0_R$(R).tsv")
    ell = d[:,1]
    luc00 = d[:,2]
    luc02 = d[:,3]
    luc20 = d[:,4]
    luc22 = d[:,5]

    w00, w02, w20, w22 = calc_wlrr_χ2303(R, ell)

    @test all(isapprox.(w00, luc00, atol=5e-11))
    @test all(isapprox.(w02, luc02, atol=5e-11))
    @test all(isapprox.(w20, luc20, atol=5e-11))
    @test w22[1] ≈ luc22[1] atol=5e-10
    @test all(isapprox.(w22[2:end], luc22[2:end], atol=5e-11))
end

function test_wl()
    test_wl_χ2303_R(1.0)
    test_wl_χ2303_R(1.1)
    test_wl_χ2303_R(0.9)
end

end # module


TestTwoFAST.test_xi()
TestTwoFAST.test_wl()


# vim: set sw=4 et sts=4 :

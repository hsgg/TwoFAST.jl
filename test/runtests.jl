#!/usr/bin/env julia


module TestTwoFAST

include("TwoFASTTestLib.jl")

using TwoFAST
using Test
using Dierckx
using .TwoFASTTestLib
using DelimitedFiles

using TwoFAST.PerformanceStats

function test_performancestats()
    t = @timed Nothing
    @show t typeof(t) TwoFAST.PerformanceStats.TimedType
    t += @timed Nothing
    t += t
    @show t
    println("t: ", t)
end

test_performancestats()


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


##################### test wljj ##############################

function test_wl_χ2303_R(R, atol)
    d = readdlm("data/wljj_chi2303.0_R$(R).tsv")
    ell = Array{Int}(d[:,1])
    luc00 = d[:,2]
    luc02 = d[:,3]
    luc20 = d[:,4]
    luc22 = d[:,5]

    w00, w02, w20, w22 = calc_wlrr_χ2303(R, ell)

    print("Testing R=$R ... ")
    @test all(isapprox.(w00, luc00, atol=atol))
    @test all(isapprox.(w02, luc02, atol=atol))
    @test all(isapprox.(w20, luc20, atol=atol))
    @test w22[1] ≈ luc22[1] atol=5e-10
    @test all(isapprox.(w22[2:end], luc22[2:end], atol=atol))
    println("passed")
end

function test_wl_ℓRR(ℓ)
    RR = [1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

    # get Lucas data
    luc = Dict()
    lucχ = readdlm("data/wljj_R0.9_ell$(ℓ).tsv")[:,1]
    for R in RR
        luc[R] = readdlm("data/wljj_R$(R)_ell$(ℓ).tsv")[:,2:end]
    end

    # get 2-FAST data
    nfast = calc_wlrr_ℓR(ℓ, RR)
    χ = nfast[1]
    wjj = nfast[2:end]
    idxmin = round(Int, length(χ) * 3/8) + 1
    idxmax = length(χ)

    # select
    s1 = (10 .<= χ .<= 1.5e4)
    s2 = (10 .<= lucχ .<= 1.5e4)

    # iterate over jj
    for jjidx=1:length(wjj)
        # iterate over R
        for i=1:length(RR)
            R = RR[i]
            w = wjj[jjidx][s1,i]
            lucw = luc[R][s2,jjidx]
            ## some debugging help:
            #println("Testing jjidx=$jjidx, R=$R...")
            #for k=1:sum(s1)
            #    @show χ[s1][k], w[k]
            #    @show lucχ[s2][k], lucw[k]
            #    diff = w[k] - lucw[k]
            #    rdiff = (w[k] - lucw[k]) / lucw[k]
            #    @show diff, rdiff
            #    @test isapprox(w[k], lucw[k], atol=5e-10, rtol=2e-4)
            #end
            print("Testing jjidx=$jjidx, R=$R...")
            @test all(isapprox.(w, lucw, atol=5e-10, rtol=2e-4))
            println(" passed")
        end
    end
end


########################## test cache ######################
function test_cache()
    # calc initial Ml cache at high ell
    RR = collect(0.6:0.1:1.1)
    cache = F21EllCache(1200, RR)
    write("out/f21ellcache", cache)
    cache2 = F21EllCache("out/f21ellcache")

    # calc full MlCache
    ell = 10:10:1200
    mlcache = MlCache(ell, "out/F21EllCache", "out/MlCache")
    write("out/MlCache", mlcache)
    mlcache = MlCache("out/MlCache")

    # ... at least a compile test ...
end


end # module


TestTwoFAST.test_xiln(0)
TestTwoFAST.test_xiln(1)
TestTwoFAST.test_xiln(2)
TestTwoFAST.test_xiln(3)
TestTwoFAST.test_xiln(4)
TestTwoFAST.test_xi_derivs()

TestTwoFAST.test_wl_χ2303_R(1.1, 1e-11)
TestTwoFAST.test_wl_χ2303_R(1.0, 2e-10)
TestTwoFAST.test_wl_χ2303_R(0.9, 1e-11)
TestTwoFAST.test_wl_χ2303_R(0.8, 1e-12)
TestTwoFAST.test_wl_χ2303_R(0.7, 1e-11)
TestTwoFAST.test_wl_χ2303_R(0.6, 1e-11)
TestTwoFAST.test_wl_χ2303_R(0.5, 2e-11)
TestTwoFAST.test_wl_χ2303_R(0.4, 1e-11)
TestTwoFAST.test_wl_χ2303_R(0.3, 1e-11)
TestTwoFAST.test_wl_χ2303_R(0.2, 1e-11)
TestTwoFAST.test_wl_χ2303_R(0.1, 1e-11)
TestTwoFAST.test_wl_ℓRR(42)
TestTwoFAST.test_wl_ℓRR(100)

TestTwoFAST.test_cache()

include("../Julia/twofast_example.jl")

# vim: set sw=4 et sts=4 :

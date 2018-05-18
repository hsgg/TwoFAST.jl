#!/usr/bin/env julia

include("TwoFASTTestLib.jl")


module TestTwoFAST

using TwoFAST
using Dierckx
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end
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

function wlrr(pk; prefix="out")
    N = 4096
    chi0 = 1e-3
    kmin = 1e-5
    kmax = 1e3
    q = 1.1
    ell = [42]  # only ell=42 for this run
    RR = [0.6, 0.7, 0.8, 0.9, 1.0]

    # calculate M_ll at high ell, result gets saved to a file:
    make_fell_lmax_cache(RR, maximum(ell), "$prefix/fell_lmax_v23.fits";
                         N=N, q=q, G=log(kmax / kmin), k0=kmin, r0=chi0)

    # calculate all M_ll, result gets saved to a file:
    tt = calcMljj(RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q,
                  fell_lmax_file="$prefix/fell_lmax_v23.fits",
                  outfile="$prefix/Ml21-cache.bin")

    # calculate wljj:
    w00 = Array{Float64}(N, length(RR))
    w02 = Array{Float64}(N, length(RR))
    function outfunc(wjj, ell, rr, RR)
        if ell == 42
            w00[:,:] = wjj[1]
            w02[:,:] = wjj[2]
        end
    end
    rr = calcwljj(pk, RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q,
                  cachefile="$prefix/Ml21-cache.bin", outfunc=outfunc)

    # store result
    writedlm("$prefix/wl42_jj00.tsv", [rr w00])
    writedlm("$prefix/wl42_jj02.tsv", [rr w02])
    println("Created '$prefix/wl42_jj00.tsv'.")
    println("Created '$prefix/wl42_jj02.tsv'.")
end


######################### call all tests ########################3

function compare_tsv(f1, f2; atol=0.0, rtol=sqrt(eps()))
    x0 = readdlm(f1)
    x1 = readdlm(f2)
    return all(isapprox.(x0, x1, atol=atol, rtol=rtol))
end

function test_all()
    mkpath("out")
    d = readdlm("data/planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat")
    pk = Spline1D(d[:,1], d[:,2])

    wlrr(pk)
    fell0, lmax0, mm0, RR0, ellmax0 = TwoFAST.read_fell_lmax("out/fell_lmax_v23.fits")
    fell1, lmax1, mm1, RR1, ellmax1 = TwoFAST.read_fell_lmax("data/fell_lmax_v23.fits")
    @test all(isapprox.(fell0, fell1))
    @test all(lmax0 .== lmax1)
    @test all(isapprox.(mm0, mm1, atol=0.0, rtol=eps(1.0)))
    #run(`cmp fell_lmax_v23.fits data/fell_lmax_v23.fits`)
    #run(`cmp Ml21-cache.bin data/Ml21-cache.bin`)
    @test compare_tsv("out/wl42_jj00.tsv", "data/wl42_jj00.tsv", atol=1e-12, rtol=1e-10)
    @test compare_tsv("out/wl42_jj02.tsv", "data/wl42_jj02.tsv", atol=1e-12, rtol=1e-10)
end

end

TestTwoFAST.test_xi()
#test_all()


# vim: set sw=4 et sts=4 :

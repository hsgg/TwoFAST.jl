#!/usr/bin/env julia


module TwoFASTTestLib

export get_quadosc_xi
export calc_2fast_xi
export calc_xiderivs
export wlrr
export calc_wlrr_χ2303

using TwoFAST
using PkSpectra


##################### ξ(r) tests ##########################

function get_quadosc_xi(ℓ)
    xi = readdlm("data/xiquadosc_plus_ell$ℓ.tsv")
    r = xi[:,1]
    xi = xi[:,2:end]
    return r, xi
end


function calc_2fast_xi(ℓ, νν)
    pk = PkSpectrum()
    kwargs = Dict(:N => 1024, :kmin => 1e-6, :kmax => 1e2, :r0 => 1e-2)
    r = xicalc(pk, 0, 0; kwargs...)[1]
    xi = [xicalc(pk, ℓ, ν; kwargs...)[2] for ν=νν]
    return r, xi
end


##################### wℓⱼⱼ′(χ,R) tests ##########################

function wlrr(RR=[0.6, 0.7, 0.8, 0.9, 1.0], ell=[42]; prefix="out", N=4096,
              chi0=1e-3, kmin=1e-5, kmax=1e3)
    mkpath(prefix)
    pk = PkSpectrum()
    q = 1.1

    # calculate M_ll at high ell, result gets saved to a file:
    make_fell_lmax_cache(RR, maximum(ell), "$prefix/fell_lmax_v23.fits";
                         N=N, q=q, G=log(kmax / kmin), k0=kmin, r0=chi0)

    # calculate all M_ll, result gets saved to a file:
    tt = calcMljj(RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q,
                  fell_lmax_file="$prefix/fell_lmax_v23.fits",
                  outfile="$prefix/Ml21-cache.bin")

    # calculate wljj:
    w00 = Array{Float64}(N, length(RR), length(ell))
    w02 = Array{Float64}(N, length(RR), length(ell))
    w20 = Array{Float64}(N, length(RR), length(ell))
    w22 = Array{Float64}(N, length(RR), length(ell))
    il = length(ell)
    function outfunc(wjj, ell, rr, RR)
        w00[:,:,il] = wjj[1]
        w02[:,:,il] = wjj[2]
        w20[:,:,il] = wjj[3]
        w22[:,:,il] = wjj[4]
        il -= 1
    end
    rr = calcwljj(pk, RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q,
                  cachefile="$prefix/Ml21-cache.bin", outfunc=outfunc)

    return w00, w02, w20, w22
end


function calc_wlrr_χ2303(R, ell)
    χ = 2303.0
    N = 1600
    kmin = 1e-5
    kmax = 1e5
    G = log(kmax / kmin)
    chi0 = χ / exp(G / 2)
    n = div(N,2) + 1
    w00, w02, w20, w22 = wlrr([R], ell; N=N, chi0=chi0, kmin=kmin, kmax=kmax)
    w00 = w00[n,1,:]
    w02 = w02[n,1,:]
    w20 = w20[n,1,:]
    w22 = w22[n,1,:]
    return w00, w02, w20, w22
end


end # module


# vim: set sw=4 et sts=4 :

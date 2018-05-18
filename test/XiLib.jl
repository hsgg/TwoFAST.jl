#!/usr/bin/env julia


module XiLib

export get_quadosc_xi
export calc_2fast_xi

using TwoFAST
using Dierckx


function get_quadosc_xi(ℓ)
    xi = readdlm("data/xiquadosc_plus_ell$ℓ.tsv")
    r = xi[:,1]
    xi = xi[:,2:end]
    return r, xi
end


function calc_2fast_xi(ℓ, νν)
    d = readdlm("data/planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat")
    pk = Spline1D(d[:,1], d[:,2])
    kwargs = Dict(:N => 1024, :kmin => 1e-6, :kmax => 1e2, :r0 => 1e-2)
    r = xicalc(pk, 0, 0; kwargs...)[1]
    xi = [xicalc(pk, ℓ, ν; kwargs...)[2] for ν=νν]
    return r, xi
end


end


# vim: set sw=4 et sts=4 :

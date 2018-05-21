#!/usr/bin/env julia


module TwoFASTTestLib

export get_quadosc_xi
export calc_2fast_xi
export calc_xiderivs

using TwoFAST
using PkSpectra


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


end


# vim: set sw=4 et sts=4 :

#!/usr/bin/env julia


module GraphicTests

using TwoFAST
using PyPlot
using Dierckx


function main()
    d = readdlm("data/planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat")
    pk = Spline1D(d[:,1], d[:,2])

    ℓ = 0

    # get 2-FAST results
    kwargs = Dict(:N => 1024, :kmin => 1e-5, :kmax => 1e3, :r0 => 1e-3)
    r0 = xicalc(pk, 0, 0; kwargs...)[1]
    xi0 = [xicalc(pk, ℓ, ν; kwargs...)[2] for ν=-2:3+ℓ]
    @show typeof(xi0)

    # get quadosc results
    xi1 = readdlm("data/xiquadosc_plus_ell$ℓ.tsv")
    r1 = xi1[:,1]

    # plot
    figure()
    for i=1:size(xi1,2)
        plot(r1, xi1[:,i], "0.5")
    end
    for i=eachindex(xi0)
        ν = (-2:3+ℓ)[i]
        plot(r0, xi0[i], "--", label="\$\\xi_$ℓ^{$ν}(r)\$")
    end
    xlim(0, 300)
    ylim(-0.01, 0.05)
    xlabel(L"r")
    ylabel(L"\xi_\ell^\nu(r)")
    legend()
end


end

GraphicTests.main()


# vim: set sw=4 et sts=4 :

#!/usr/bin/env julia

include("XiLib.jl")

module GraphicTests

using TwoFAST
using XiLib
using PyCall
using PyPlot
using Dierckx
@pyimport matplotlib.gridspec as gridspec


function plot_xi_rdiff(ℓ)
    νν = -2:3
    r0, xi0 = get_quadosc_xi(ℓ)
    r1, xi1 = calc_2fast_xi(ℓ, νν)

    # plot
    figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
    gs[:update](hspace=0)
    ax1 = subplot(gs[1])
    ax2 = subplot(gs[2], sharex=ax1)
    setp(ax1[:get_xticklabels](), visible=false)
    ax1[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)
    ax2[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)
    for i=eachindex(xi1)
        ν = νν[i]
        s0 = Spline1D(r0, xi0[:,i])
        s1 = Spline1D(r1, xi1[i])
        r = r0[1:10:end]
        diff = s1(r) - s0(r)
        rdiff = diff ./ s0(r)
        ax1[:plot](r0, xi0[:,i], "0.5", alpha=0.75, zorder=-1)
        ax1[:plot](r1, xi1[i], "--", label="\$\\xi_$ℓ^{$ν}(r)\$")
        ax2[:plot](r, rdiff, "--")
    end
    ax1[:set_xlim](0, 200)
    ax1[:set_ylim](-0.01, 0.05)
    ax2[:set_ylim](-0.00019, 0.00019)
    ax2[:set_xlabel](L"r")
    ax1[:set_ylabel](L"\xi_\ell^\nu(r)")
    ax2[:set_ylabel](L"\Delta\xi/\xi")
    ax1[:legend]()
    tight_layout()
end


end

close("all")
GraphicTests.plot_xi_rdiff(0)
GraphicTests.plot_xi_rdiff(1)
GraphicTests.plot_xi_rdiff(2)
GraphicTests.plot_xi_rdiff(3)
GraphicTests.plot_xi_rdiff(4)


# vim: set sw=4 et sts=4 :

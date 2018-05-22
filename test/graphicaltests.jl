#!/usr/bin/env julia

include("PkSpectra.jl")
include("TwoFASTTestLib.jl")

module GraphicTests

using TwoFAST
using TwoFASTTestLib
using PyCall
using PyPlot
using Dierckx
@pyimport matplotlib.gridspec as gridspec


##################### ξ(r) tests ##########################

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


function plot_xi_derivs()
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

    r = rq[1:10:end]

    # plot
    figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
    gs[:update](hspace=0)
    ax1 = subplot(gs[1])
    ax2 = subplot(gs[2], sharex=ax1)
    setp(ax1[:get_xticklabels](), visible=false)
    ax1[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)
    ax2[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)

    # ξ
    diff = ξ(r) - ξq(r)
    rdiff = diff ./ ξq(r)
    ax1[:plot](rq, ξq(rq), "0.5", alpha=0.75, zorder=-1)
    ax1[:plot](r0, ξ(r0), "--", label="\$\\xi(r)\$")
    ax2[:plot](r, rdiff, "--")

    # ξ′
    diff = ξ′(r) - ξq′(r)
    rdiff = diff ./ ξq′(r)
    ax1[:plot](rq, ξq′(rq), "0.5", alpha=0.75, zorder=-1)
    ax1[:plot](r0, ξ′(r0), "--", label="\$\\xi'(r)\$")
    ax2[:plot](r, rdiff, "--")

    # ξ″
    diff = ξ″(r) - ξq″(r)
    rdiff = diff ./ ξq″(r)
    ax1[:plot](rq, ξq″(rq), "0.5", alpha=0.75, zorder=-1)
    ax1[:plot](r0, ξ″(r0), "--", label="\$\\xi''(r)\$")
    ax2[:plot](r, rdiff, "--")

    ax1[:set_xlim](0, 200)
    ax1[:set_ylim](-0.001, 0.003)
    ax2[:set_ylim](-0.00049, 0.00049)
    ax2[:set_xlabel](L"r")
    ax1[:set_ylabel](L"\xi(r)")
    ax2[:set_ylabel](L"\Delta\xi/\xi")
    ax1[:legend]()
    tight_layout()
end


##################### wℓⱼⱼ′(χ,R) tests ##########################

function plot_cl(R=1.0, xmax=1200, y2max=1.5e-6)
    println("plot_cl(R=$R)...")

    d = readdlm("data/wljj_chi2303.0_R$(R).tsv")
    lucell = d[:,1]
    luc00 = d[:,2]
    luc02 = d[:,3]
    luc20 = d[:,4]
    luc22 = d[:,5]

    ell = lucell
    w00, w02, w20, w22 = calc_wlrr_χ2303(R, ell)

    fac = 1
    #fac = ell
    #fac = ell * (ell + 1) / 2 / np.pi

    figure()
    gs = gridspec.GridSpec(3, 1, height_ratios=[4,1,1])
    gs[:update](hspace=0)
    ax1 = subplot(gs[1])
    ax2 = subplot(gs[2], sharex=ax1)
    ax3 = subplot(gs[3], sharex=ax1)
    setp(ax1[:get_xticklabels](), visible=false)
    ax1[:plot](lucell, fac * luc00, color="0.55")
    ax1[:plot](lucell, fac * luc02, color="0.65")
    ax1[:plot](lucell, fac * luc20, color="0.75")
    ax1[:plot](lucell, fac * luc22, color="0.85")
    ax1[:plot](ell, fac * w00, label=L"(j,j')=(0,0)", ls="--")
    ax1[:plot](ell, fac * w02, label=L"(j,j')=(0,2)", ls="--")
    ax1[:plot](ell, fac * w20, label=L"(j,j')=(2,0)", ls="--")
    ax1[:plot](ell, fac * w22, label=L"(j,j')=(2,2)", ls="--")
    ax1[:hlines](0.0, 0, 1200, color="0.85")
    ax1[:set_ylabel](L"$w_{\ell,jj'}$")
    ax1[:set_xlabel](L"$\ell$")
    ax1[:ticklabel_format](styl="sci", scilimits=(-1,4))
    ax1[:legend](loc="upper right", bbox_to_anchor=(1.02, 1.04))

    ax2[:plot](ell, w00 - luc00, label=L"diff, jj'=00")
    ax2[:plot](ell, w02 - luc02, label=L"diff, jj'=02")
    ax2[:plot](ell, w20 - luc20, label=L"diff, jj'=20")
    ax2[:plot](ell, w22 - luc22, label=L"diff, jj'=22")
    ax2[:set_ylabel](L"$\Delta w$")
    ax2[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)

    ax3[:plot](ell, (w00 - luc00) ./ luc00, label="rdiff, jj'=00", ls="--")
    ax3[:plot](ell, (w02 - luc02) ./ luc02, label="rdiff, jj'=02", ls="--")
    ax3[:plot](ell, (w20 - luc20) ./ luc20, label="rdiff, jj'=20", ls="--")
    ax3[:plot](ell, (w22 - luc22) ./ luc22, label="rdiff, jj'=22", ls="--")
    ax3[:set_ylabel](L"$\Delta w / w$")
    ax3[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)

    suptitle(L"$R="*string(R)*L"$", x=0.55)
    tight_layout()
    #savefig("wl22_R$R.pdf")
end



end # module

close("all")
#GraphicTests.plot_xi_rdiff(0)
#GraphicTests.plot_xi_rdiff(1)
#GraphicTests.plot_xi_rdiff(2)
#GraphicTests.plot_xi_rdiff(3)
#GraphicTests.plot_xi_rdiff(4)
#GraphicTests.plot_xi_derivs()

GraphicTests.plot_cl(1.0)
GraphicTests.plot_cl(1.1)
GraphicTests.plot_cl(0.9)


# vim: set sw=4 et sts=4 :

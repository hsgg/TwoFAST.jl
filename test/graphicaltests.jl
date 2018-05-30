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


########################### utilities #######################

function logplot(x, y; ax=gca(), kwargs...)
    select = isfinite.(y)
    x = x[select]
    y = y[select]

    # plot positive values
    l = ax[:plot](x, y; kwargs...)

    kwargs = Dict{Symbol,Any}(kwargs)
    kwargs[:color] = l[1][:get_color]()
    kwargs[:alpha] = l[1][:get_alpha]()
    kwargs[:linestyle] = l[1][:get_linestyle]()
    try delete!(kwargs, :label) end
    try delete!(kwargs, :ls) end


    if kwargs[:linestyle] == "-"
        kwargs[:linestyle] = "--"
    elseif kwargs[:linestyle] == "--"
        kwargs[:linestyle] = ":"
    end

    ax[:plot](x, -y; kwargs...)
end

##################### ξ(r) tests ##########################

function plot_xi_rdiff(ℓ)
    νν = -2:3
    r0, xi0 = get_quadosc_xi(ℓ)
    r1, xi1 = calc_2fast_xi(ℓ, νν)

    # plot
    figure()
    gs = gridspec.GridSpec(3, 1, height_ratios=[3,1,1])
    gs[:update](hspace=0)
    ax1 = subplot(gs[1])
    ax2 = subplot(gs[2], sharex=ax1)
    ax3 = subplot(gs[3], sharex=ax1)
    setp(ax1[:get_xticklabels](), visible=false)
    setp(ax2[:get_xticklabels](), visible=false)
    ax1[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)
    ax2[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)
    ax3[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)
    for i=eachindex(xi1)
        ν = νν[i]
        s0 = Spline1D(r0, xi0[:,i])
        s1 = Spline1D(r1, xi1[i])
        r = r0[1:10:end]
        diff = s1(r) - s0(r)
        rdiff = diff ./ s0(r)
        ax1[:plot](r0, xi0[:,i], "0.5", alpha=0.75, zorder=-1)
        ax1[:plot](r1, xi1[i], "--", label="\$\\xi_$ℓ^{$ν}(r)\$")
        ax2[:plot](r, diff, "--")
        ax3[:plot](r, rdiff, "--")
    end
    ax1[:set_xlim](0, 200)
    ax1[:set_ylim](-0.01, 0.05)
    ax1[:set_ylabel](L"\xi_\ell^\nu(r)")
    ax1[:legend]()
    ax2[:set_ylabel](L"\Delta\xi")
    ax2[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)
    #ax3[:set_ylim](-0.00019, 0.00019)
    ax3[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)
    ax3[:set_xlabel](L"r")
    ax3[:set_ylabel](L"\Delta\xi/\xi")
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
    gs = gridspec.GridSpec(3, 1, height_ratios=[3,1,1])
    gs[:update](hspace=0)
    ax1 = subplot(gs[1])
    ax2 = subplot(gs[2], sharex=ax1)
    ax3 = subplot(gs[3], sharex=ax1)
    setp(ax1[:get_xticklabels](), visible=false)
    setp(ax2[:get_xticklabels](), visible=false)
    ax1[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)
    ax2[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)
    ax3[:hlines](0.0, 0.0, 200, color="0.85", zorder=-2)

    # ξ
    diff = ξ(r) - ξq(r)
    rdiff = diff ./ ξq(r)
    ax1[:plot](rq, ξq(rq), "0.5", alpha=0.75, zorder=-1)
    ax1[:plot](r0, ξ(r0), "--", label="\$\\xi(r)\$")
    ax2[:plot](r, diff, "--")
    ax3[:plot](r, rdiff, "--")

    # ξ′
    diff = ξ′(r) - ξq′(r)
    rdiff = diff ./ ξq′(r)
    ax1[:plot](rq, ξq′(rq), "0.5", alpha=0.75, zorder=-1)
    ax1[:plot](r0, ξ′(r0), "--", label="\$\\xi'(r)\$")
    ax2[:plot](r, diff, "--")
    ax3[:plot](r, rdiff, "--")

    # ξ″
    diff = ξ″(r) - ξq″(r)
    rdiff = diff ./ ξq″(r)
    ax1[:plot](rq, ξq″(rq), "0.5", alpha=0.75, zorder=-1)
    ax1[:plot](r0, ξ″(r0), "--", label="\$\\xi''(r)\$")
    ax2[:plot](r, diff, "--")
    ax3[:plot](r, rdiff, "--")

    ax1[:set_xlim](0, 200)
    ax1[:set_ylim](-0.001, 0.003)
    ax2[:set_ylim](-0.00049, 0.00049)
    ax2[:set_xlabel](L"r")
    ax2[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)
    ax1[:set_ylabel](L"\xi(r)")
    ax2[:set_ylabel](L"\Delta\xi")
    ax3[:set_ylabel](L"\Delta\xi/\xi")
    ax3[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)
    ax1[:legend]()
    tight_layout()
end


##################### wℓⱼⱼ′(χ,R) tests ##########################
# along ℓ
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
    setp(ax2[:get_xticklabels](), visible=false)
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
    ax1[:ticklabel_format](styl="sci", scilimits=(-1,4))
    ax1[:legend](loc="upper right")

    ax2[:plot](ell, w00 - luc00, label=L"diff, jj'=00")
    ax2[:plot](ell, w02 - luc02, label=L"diff, jj'=02")
    ax2[:plot](ell, w20 - luc20, label=L"diff, jj'=20")
    ax2[:plot](ell, w22 - luc22, label=L"diff, jj'=22")
    ax2[:set_ylabel](L"$\Delta w$")
    ax2[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)
    #ax2[:get_yaxis]()[:get_offset_text]()[:set_y](0.5)
    #ax2[:yaxis]()[:offsetText]()

    ax3[:plot](ell, (w00 - luc00) ./ luc00, label="rdiff, jj'=00", ls="--")
    ax3[:plot](ell, (w02 - luc02) ./ luc02, label="rdiff, jj'=02", ls="--")
    ax3[:plot](ell, (w20 - luc20) ./ luc20, label="rdiff, jj'=20", ls="--")
    ax3[:plot](ell, (w22 - luc22) ./ luc22, label="rdiff, jj'=22", ls="--")
    ax3[:set_ylabel](L"$\Delta w / w$")
    ax3[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)
    ax3[:set_xlabel](L"$\ell$")

    suptitle(L"$R="*string(R)*L"$", x=0.55)
    tight_layout()
    #savefig("wl22_R$R.pdf")
end

# along χ
function plot_cl_χ(ℓ=42, jjidx=1, y2max=1.5e-6)
    println("plot_cl_χ(ℓ=$ℓ, jjidx=$jjidx)...")
    RR = [1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    #RR = [0.9]

    # get Lucas data
    luc = Dict()
    lucχ = readdlm("data/wljj_R0.9_ell$(ℓ).tsv")[:,1]
    for R in RR
        luc[R] = readdlm("data/wljj_R$(R)_ell$(ℓ).tsv")[:,jjidx+1]
    end

    # get 2-FAST data
    tfast = calc_wlrr_ℓR(ℓ, RR)
    χ = tfast[1]
    wjj = tfast[jjidx+1]

    fac = 1
    #fac = ℓ
    #fac = ℓ * (ℓ + 1) / 2π

    figure(figsize=[6.4, 4.8]*0.73/0.49)
    gs = gridspec.GridSpec(3, 1, height_ratios=[4,1,1])
    gs[:update](hspace=0)
    ax1 = subplot(gs[1])
    ax2 = subplot(gs[2], sharex=ax1)
    ax3 = subplot(gs[3], sharex=ax1)
    setp(ax1[:get_xticklabels](), visible=false)
    setp(ax2[:get_xticklabels](), visible=false)
    for i=1:length(RR)
        logplot(lucχ, fac * luc[RR[i]], ax=ax1, color="0.55")
        logplot(χ, fac * wjj[:,i], ax=ax1, label="\$R=$(RR[i])\$", ls="--")
    end
    ax1[:set_ylabel](L"$w_{\ell,jj'}$")
    ax1[:ticklabel_format](styl="sci", scilimits=(-1,4))
    ax1[:legend](loc="upper right")
    ax1[:set_yscale]("log")

    for i=1:length(RR)
        w = Spline1D(χ, wjj[:,i])(lucχ)
        ax2[:plot](lucχ, w - luc[RR[i]], label="diff, \$R=$(RR[i])\$", ls="--")
        ax3[:plot](lucχ, (w - luc[RR[i]]) ./ luc[RR[i]], label="rdiff, \$R=$(RR[i])\$", ls="--")
    end
    ax2[:set_ylabel](L"$\Delta w$")
    ax2[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)
    ax3[:set_ylabel](L"$\Delta w / w$")
    ax3[:ticklabel_format](style="sci", axis="y", scilimits=(0,0), useOffset=false)
    ax3[:set_xlabel](L"$\chi$")
    xscale("log")
    xlim(1e1, 1.5e4)

    suptitle("\$\\ell=$ℓ\$, \$(j,j')=($jjidx)\$", y=1.01)
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
GraphicTests.plot_cl(0.8)
GraphicTests.plot_cl(0.7)
GraphicTests.plot_cl(0.6)
GraphicTests.plot_cl(0.5)
GraphicTests.plot_cl(0.4)
GraphicTests.plot_cl(0.3)
GraphicTests.plot_cl(0.2)
GraphicTests.plot_cl(0.1)

GraphicTests.plot_cl_χ(42, 1)
GraphicTests.plot_cl_χ(42, 2)
GraphicTests.plot_cl_χ(42, 3)
GraphicTests.plot_cl_χ(42, 4)

show()

# vim: set sw=4 et sts=4 :

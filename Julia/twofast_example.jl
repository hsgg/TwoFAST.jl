#!/usr/bin/env julia


using TwoFAST

# for Spline1D (add it with Pkg.add("Dierckx")):
using Dierckx


# one spherical bessel:
function xiln(pk; N=1024, fname="xi.tsv")
	kmax = 1e3
	kmin=1e-5
	r0 = 1e-3

	print("xi, l=0, nu=0: ")
	@time r00, xi00 = xicalc(pk, 0, 0; N=N, kmin=kmin, kmax=kmax, r0=r0)

	print("xi, l=0, nu=-2:")
	@time r, xi0m2 = xicalc(pk, 0, -2; N=N, kmin=kmin, kmax=kmax, r0=r0)

	writedlm(fname, [r00 xi00 xi0m2])
	println("Created '$fname'.")
end


# two spherical bessel:
function wlrr(pk)
	N = 4096
	chi0 = 1e-3
	kmin = 1e-5
	kmax = 1e3
	q = 1.1
	ell = [42]  # only ell=42 for this run
	RR = [0.6, 0.7, 0.8, 0.9, 1.0]

	# calculate M_ll at high ell:
	make_fell_lmax_cache(RR, maximum(ell); N=N, q=q, G=log(kmax / kmin), k0=kmin, r0=chi0)

	# calculate all M_ll:
	tt = calcMljj(RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q)

	# store result here:
	w00 = Array{Float64}(N, length(RR))
	w02 = Array{Float64}(N, length(RR))
	function outfunc(wjj, ell, rr, RR)
		if ell == 42
			w00[:,:] = wjj[1]
			w02[:,:] = wjj[2]
		end
	end
	rr = calcwljj(pk, RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q,
		outfunc=outfunc)
	writedlm("wl42_jj00.tsv", [rr w00])
	writedlm("wl42_jj02.tsv", [rr w02])
	println("Created 'wl42_jj00.tsv'.")
	println("Created 'wl42_jj02.tsv'.")
end


function main()
	d = readdlm("../planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat")
	pk = Spline1D(d[:,1], d[:,2])

	# one spherical bessel function:
	xiln(pk)

	# two spherical bessel functions:
	wlrr(pk)
end



main()

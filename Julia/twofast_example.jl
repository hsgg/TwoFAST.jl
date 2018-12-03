#!/usr/bin/env julia


using TwoFAST

# for Spline1D (add it with Pkg.add("Dierckx")):
using Dierckx
using DelimitedFiles


# one spherical bessel:
function xiln(pk; N=1024, fname="out/xi.tsv")
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
	ell = [42, 1000]  # only ell=42 for this run
	RR = [0.6, 0.7, 0.8, 0.9, 0.998, 1.0]

	# calculate M_ll at high ell:
	f21cache = F21EllCache(maximum(ell), RR, N; q=q, kmin=kmin, kmax=kmax,
			       χ0=chi0)
	write("out/F21EllCache", f21cache)

	# calculate all M_ll:
	mlcache = MlCache(ell, "out/F21EllCache", "out/MlCache")
	write("out/MlCache", mlcache)

	# store result here:
	w00 = Array{Float64}(undef, N, length(RR))
	w02 = Array{Float64}(undef, N, length(RR))
	function outfunc(wjj, ell, rr, RR)
		if ell == 42
			w00[:,:] = wjj[1]
			w02[:,:] = wjj[2]
		end
	end
	rr = calcwljj(pk, RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q,
		outfunc=outfunc, cachefile="out/MlCache/MlCache.bin")
	writedlm("out/wl42_jj00.tsv", [rr w00])
	writedlm("out/wl42_jj02.tsv", [rr w02])
	println("Created 'out/wl42_jj00.tsv'.")
	println("Created 'out/wl42_jj02.tsv'.")
end


function main()
	d = readdlm((@__DIR__)*"/../test/data/planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat", comments=true)
	pk = Spline1D(d[:,1], d[:,2])

	# one spherical bessel function:
	xiln(pk)

	# two spherical bessel functions:
	wlrr(pk)
end



main()

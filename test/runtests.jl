#!/usr/bin/env julia

using TwoFAST
using Dierckx
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end


############### test xicalc #######################

function xiln(pk; N=1024, fname="xi.tsv")
	kmax = 1e3
	kmin = 1e-5
	r0 = 1e-3

	print("xi, l=0, nu=0: ")
	@time r00, xi00 = xicalc(pk, 0, 0; N=N, kmin=kmin, kmax=kmax, r0=r0)

	print("xi, l=0, nu=-2:")
	@time r, xi0m2 = xicalc(pk, 0, -2; N=N, kmin=kmin, kmax=kmax, r0=r0)

	writedlm(fname, [r00 xi00 xi0m2])
	println("Created '$fname'.")
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

	xiln(pk; fname="out/xi.tsv")
	@test compare_tsv("out/xi.tsv", "data/xi.tsv", atol=1e-10, rtol=1e-10)

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

test_all()

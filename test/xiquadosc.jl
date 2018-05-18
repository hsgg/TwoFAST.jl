#!/usr/bin/env julia

using Hwloc
addprocs(Hwloc.num_physical_cores())
println("hostname: ", gethostname())
println("julia version: ", VERSION)
println("nprocs: ", nprocs())
println("nworkers: ", nworkers())


@everywhere dir = "/home/hsgg/research/hgebhardt/code/mypy"
try
	@everywhere include("$dir/pkspec.jl")
	@everywhere include("$dir/bisect.jl")
	@everywhere include("$dir/quadosc.jl")
	@everywhere include("$dir/sphbes/sphbes.jl")
catch
	println("Cannot load essential libraries from `$dir`")
	println("They are not included because the license is unclear to me.")
	println("What is the numerical recipes license?")
	error("Cannot load essential libraries.")
end

@everywhere using PkSpectrum
@everywhere using QuadOsc
@everywhere using SphBes
@everywhere using QuadGK: quadgk



@everywhere function integrand(lnk, r, pwr, ell=0, nu=0)
    if lnk > 700
        return 0.0
    end
    k = exp(lnk)
    if k > 1e50
        return 0.0
    end
    res = sphbesj(k*r, ell) * k^(3-nu) * pwr(k)
    #println("lnk=$lnk\tk=$k\tnu=$nu\tres=$res")
    if !isfinite(res)
        error("not finite!")
    end
    return res
end


@everywhere function corrfunc(r::Number, pwr, ell=0, nu=0; reltol=1e-10, order=127)
	if r == 0
		return 0.0
	end
	integ(lnk) = integrand(lnk, r, pwr, ell, nu)
	pivot = 0.0
	#print("computing r=$r... ")
	I0, E0 = quadgk(integ, -Inf, pivot, reltol=reltol, order=order)
	I1, E1 = quadosc(integ, pivot, Inf, n -> log(n * pi / r), beta=1.0,
		reltol=reltol, order=order)
	factor = 1 / (2 * pi^2 * r^nu)
	return factor * (I0 + I1), factor * (E0 + E1)
end

function corrfunc(rr, pwr, ell=0, nu=0; reltol=1e-10, order=127)
	print("ξ(r), ℓ=$ell, ν=$nu:")
	if (nu == -2 && ell ∈ [0, 1, 2, 3, 4] ||
	    nu == -1 && ell ∈ [1])
		@time nothing
		return fill(NaN, length(rr)), fill(NaN, length(rr))
	end
	@time tmp = pmap(r -> corrfunc(r, pwr, ell, nu, reltol=reltol, order=order), rr)
	xi = Array{Float64}(length(rr))
	xie = Array{Float64}(length(rr))
	for i=1:length(rr)
		xi[i] = tmp[i][1]
		xie[i] = tmp[i][2]
	end
	return xi, xie
end



function main()
	fname = "data/planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat"
	pkin = PkFile(fname)#; ns=0.9672)
	r = linspace(1.0, 200.0, Int(199/0.1) + 1)

	# Note: Some of the combinations of ell and nu below should converge,
	# but take a long, long time.

	print("xi, l=0, nu=0: ")
	@time xi00, xi00e = corrfunc(r, pkin, 0, 0, reltol=1e-13, order=511)
	# Note: reltol=1e-14 doesn't improve it, but takes 12 hours instead of 2 minutes!

	#print("xi, l=0, nu=-2:")
	#@time xi0m2, xi0m2e = corrfunc(r, pkin, 0, -2)  # fails
	xi0m2 = NaN * zeros(length(r))
	xi0m2e = Inf * zeros(length(r))

	#print("xi, l=1, nu=-1:")
	#@time xi1m1, xi1m1e = corrfunc(r, pkin, 1, -1)  # takes super-long
	xi1m1 = NaN * zeros(length(r))
	xi1m1e = Inf * zeros(length(r))

	#print("xi, l=1, nu=1: ")
	#@time xi11, xi11e = corrfunc(r, pkin, 1, 1)
	xi11 = NaN * zeros(length(r))
	xi11e = Inf * zeros(length(r))

	print("xi, l=1, nu=3: ")
	@time xi13, xi13e = corrfunc(r, pkin, 1, 3)
	#xi13 = NaN * zeros(length(r))
	#xi13e = Inf * zeros(length(r))

	#print("xi, l=2, nu=-2:")
	#@time xi2m2 xi2m2e = corrfunc(r, pkin, 2, -2)  # fails
	xi2m2 = NaN * zeros(length(r))
	xi2m2e = Inf * zeros(length(r))

	print("xi, l=2, nu=0: ")
	@time xi20, xi20e = corrfunc(r, pkin, 2, 0)
	#xi20 = NaN * zeros(length(r))
	#xi20e = Inf * zeros(length(r))

	#print("xi, l=3, nu=-1:")
	#@time xi3m1, xi3m1e = corrfunc(r, pkin, 3, -1)
	xi3m1 = NaN * zeros(length(r))
	xi3m1e = Inf * zeros(length(r))

	#print("xi, l=3, nu=1: ")
	#@time xi31, xi31e = corrfunc(r, pkin, 3, 1)
	xi31 = NaN * zeros(length(r))
	xi31e = Inf * zeros(length(r))

	print("xi, l=4, nu=0: ")
	@time xi40, xi40e = corrfunc(r, pkin, 4, 0)
	#xi40 = NaN * zeros(length(r))
	#xi40e = Inf * zeros(length(r))


	writedlm("data/xiquadosc.dat", [r xi00 xi0m2 xi1m1 xi11 xi13 xi2m2 xi20 xi3m1 xi31 xi40])
	writedlm("data/xiquadoscerr.dat", [r xi00e xi0m2e xi1m1e xi11e xi13e xi2m2e xi20e xi3m1e xi31e xi40e])
end


function more_xi(ℓ=0)
	fname = "data/planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat"
	pkin = PkFile(fname)#; ns=0.9672)
	r = linspace(1.0, 200.0, Int(199/0.1) + 1)

	xiℓ = [corrfunc(r, pkin, ℓ, ν)[1] for ν=-2:3]

	writedlm("data/xiquadosc_plus_ell$ℓ.tsv", [r xiℓ...])
end


#main()
more_xi(0)
more_xi(1)
more_xi(2)
more_xi(3)
more_xi(4)

#!/usr/bin/env julia




#using Hwloc
#num_physical_cores = Hwloc.num_physical_cores()
num_physical_cores = 2

using Distributed
println("nprocs: ", nprocs())
println("nworkers: ", nworkers())
addprocs(num_physical_cores + 1 - nprocs())
println("hostname: ", gethostname())
println("julia version: ", VERSION)
println("nprocs: ", nprocs())
println("nworkers: ", nworkers())


@everywhere include("PkSpectra.jl")
@everywhere include("$(homedir())/research/hgebhardt/myjl/include.jl")

@everywhere using .PkSpectra
@everywhere using QuadOsc
@everywhere using SphBes
@everywhere using QuadGK
@everywhere using DelimitedFiles



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


@everywhere function corrfunc(r::Number, pwr, ell=0, nu=0; rtol=1e-10, order=127)
	if r == 0
		return 0.0
	end
	integ(lnk) = integrand(lnk, r, pwr, ell, nu)
	pivot = 0.0
	#print("computing r=$r... ")
	I0, E0 = quadgk(integ, -Inf, pivot, rtol=rtol, order=order)
	I1, E1 = quadosc(integ, pivot, Inf, n -> log(n * pi / r), beta=1.0,
		rtol=rtol, order=order)
	factor = 1 / (2 * pi^2 * r^nu)
	return factor * (I0 + I1), factor * (E0 + E1)
end

function corrfunc(rr, pwr, ell=0, nu=0; rtol=1e-10, order=127)
	print("ξ(r), ℓ=$ell, ν=$nu: ")
	if (nu == -2 && ell ∈ [0, 1, 2, 3, 4] ||
	    nu == -1 && ell ∈ [1])
		println("Skipping...")
		return fill(NaN, length(rr)), fill(NaN, length(rr))
	end
	@time tmp = pmap(r -> corrfunc(r, pwr, ell, nu, rtol=rtol, order=order), rr)
	xi = Array{Float64}(undef, length(rr))
	xie = Array{Float64}(undef, length(rr))
	for i=1:length(rr)
		xi[i] = tmp[i][1]
		xie[i] = tmp[i][2]
	end
	return xi, xie
end


function more_xi(ℓ=0)
	pkin = PkSpectrum()
	r = range(1.0, stop=200.0, length=Int(199/0.1) + 1)

	xiℓ = [corrfunc(r, pkin, ℓ, ν)[1] for ν=-2:3]

	writedlm("data/xiquadosc_plus_ell$ℓ.tsv", [r xiℓ...])
end


more_xi(0)
more_xi(1)
more_xi(2)
more_xi(3)
more_xi(4)

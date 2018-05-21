#!/usr/bin/env julia
##!/home/hsgg/hetdex/julia

# Written by Henry Gebhardt (2016-2017)

using Hwloc
#addprocs(Hwloc.num_physical_cores())
addprocs(4)
println("hostname: ", gethostname())
println("julia version: ", VERSION)
println("nprocs: ", nprocs())
println("nworkers: ", nworkers())


@everywhere include("PkSpectra.jl")

@everywhere dir = "$(homedir())/research/hgebhardt/code/mypy"
try
	@everywhere include("$dir/bisect.jl")
	@everywhere include("$dir/sphbes/sphbes.jl")
	@everywhere include("$dir/quad_jar_jbt.jl")
catch
	println("Cannot load essential libraries from `$dir`")
	println("They are not included because the license is unclear to me.")
	println("What is the numerical recipes license?")
	error("Cannot load essential libraries.")
end


#=
# trapezoidal
function wllrr(k3pk, dlnk, kr1, kr2, ell1, ell2)
	@inline kern(i) = k3pk[i] * sphbesj(kr1[i], ell1) * sphbesj(kr2[i], ell2)
	s = kern(1) / 2
	for i in 2:length(k3pk)-1
		s += kern(i)
	end
	s += kern(length(k3pk)) / 2
	s * dlnk * (2/pi)
end
function wllrr(ell1, ell2, r1::Number, r2::Number, pk, beta=3)
	N = 2^20
	k = logspace(-3, 3, N)
	dlnk = (log(maximum(k)) - log(minimum(k))) / N
	k3pk = k.^beta .* pk(k)
	kr1 = k .* r1
	kr2 = k .* r2
	I = wllrr(k3pk, dlnk, kr1, kr2, ell1, ell2)
	return I, 0.0
end
=#


@everywhere module WlLucas

using PkSpectra
using Quadjarjbt

# Lucas 1995
function wllrr(ell1::Integer, ell2::Integer, r1::Number, r2::Number, pwr, beta=3; reltol=1e-10)
	function f(lnk)
		if lnk > 236.0  # = log(cbrt(FLOAT_MAX))
			ret = 0.0
		else
			k = exp(lnk)
			ret = k^beta * pwr(k)
		end
		if !isfinite(ret)
			println("lnk=$lnk, ret=$ret")
		end
		return ret
	end
	t = @elapsed I, E = quad_jar_jbt_log(f, -Inf, Inf, ell1, ell2, r1, r2;
		n1=ell1>3?ell1^2:3,
		n2=ell2>3?ell2^2:3,
		reltol=reltol, abstol=1e-20, order=511)
	I *= 2/pi
	E *= 2/pi  # Note: Error E is not reliable
	println("wl(ℓ₁=$ell1, ℓ₂=$ell2, r₁=$r1, r₂=$r2) = $I: took $t sec")
	I
end


# calc along ℓ, Δℓ=-4,-2,0,2,4
function calc_wlχR(χ, R, elllist=[2:100,112,125,150,200,300,400,500,600,700,800,900,1000,1200])
	pk = PkSpectrum()

	ell = Int[]
	wlm4 = Float64[]
	wlm2 = Float64[]
	wl0  = Float64[]
	wlp2 = Float64[]
	wlp4 = Float64[]
	ellgood = Bool[]
	for ells ∈ elllist
		ellrange = minimum(ells)-2:maximum(ells)+2
		append!(ell, ellrange)
		append!(wlm4, pmap(ℓ->wllrr(ℓ, ℓ-4, χ, R*χ, pk), ellrange))
		append!(wlm2, pmap(ℓ->wllrr(ℓ, ℓ-2, χ, R*χ, pk), ellrange))
		append!(wl0 , pmap(ℓ->wllrr(ℓ, ℓ,   χ, R*χ, pk), ellrange))
		append!(wlp2, pmap(ℓ->wllrr(ℓ, ℓ+2, χ, R*χ, pk), ellrange))
		append!(wlp4, pmap(ℓ->wllrr(ℓ, ℓ+4, χ, R*χ, pk), ellrange))
		append!(ellgood, [ℓ ∈ ells for ℓ ∈ ellrange])
	end

	fname = "data/wl_chi$(χ)_R$(R).tsv"
	writedlm(fname, [ell wlm4 wlm2 wl0 wlp2 wlp4 ellgood])
	println("Created '$fname'.")
end

end # module
using WlLucas


WlLucas.calc_wlχR(2303.0, 1.0, [2:1200])

#χ = logspace(0, 4)

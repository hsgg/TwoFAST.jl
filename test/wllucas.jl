#!/usr/bin/env julia
##!/home/hsgg/hetdex/julia

# Written by Henry Gebhardt (2016-2017)

using Hwloc
println("nprocs: ", nprocs())
println("nworkers: ", nworkers())
addprocs(Hwloc.num_physical_cores() + 1 - nprocs())
#addprocs(4 + 1 - nprocs())
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
function wllrr(ell1::Integer, ell2::Integer, r1::Number, r2::Number, pwr, n=0; reltol=1e-10)
	beta = 3 + n
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
		reltol=reltol, abstol=1e-14, order=511)
	I *= 2/pi
	E *= 2/pi  # Note: Error E is not reliable
	println("wl(ℓ₁=$ell1, ℓ₂=$ell2, r₁=$r1, r₂=$r2) = $I: took $t sec")
	I
end


# calc along ℓ, Δℓ=0,±2,±4
function calc_wldlχR(χ, R, elllist=[2:100,112,125,150,200,300,400,500,600,700,800,900,1000,1200],
		     n=0, fname="data/wldl_chi$(χ)_R$(R).tsv")
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
		append!(wlm4, pmap(ℓ->wllrr(ℓ, ℓ-4, χ, R*χ, pk, n), ellrange))
		append!(wlm2, pmap(ℓ->wllrr(ℓ, ℓ-2, χ, R*χ, pk, n), ellrange))
		append!(wl0 , pmap(ℓ->wllrr(ℓ, ℓ,   χ, R*χ, pk, n), ellrange))
		append!(wlp2, pmap(ℓ->wllrr(ℓ, ℓ+2, χ, R*χ, pk, n), ellrange))
		append!(wlp4, pmap(ℓ->wllrr(ℓ, ℓ+4, χ, R*χ, pk, n), ellrange))
		append!(ellgood, [ℓ ∈ ells for ℓ ∈ ellrange])
	end

	writedlm(fname, [ell wl0 wlm2 wlp2 wlm4 wlp4 ellgood])
	println("Created '$fname'.")
end


# calc along χ, Δℓ=0,±2,±4
function calc_wldlRℓ(R, ℓ, n, χrange; fname="data/wldl_R$(R)_ell$(ℓ).tsv")
	pk = PkSpectrum()
	wl0  = pmap(χ->wllrr(ℓ, ℓ,   χ, R*χ, pk, n), χrange)
	wlm2 = pmap(χ->wllrr(ℓ, ℓ-2, χ, R*χ, pk, n), χrange)
	wlp2 = pmap(χ->wllrr(ℓ, ℓ+2, χ, R*χ, pk, n), χrange)
	wlm4 = pmap(χ->wllrr(ℓ, ℓ-4, χ, R*χ, pk, n), χrange)
	wlp4 = pmap(χ->wllrr(ℓ, ℓ+4, χ, R*χ, pk, n), χrange)
	writedlm(fname, [χrange wl0 wlm2 wlp2 wlm4 wlp4])
	println("Created '$fname'.")
end


function wljj_dl(ell, j1, j2, wldl; abs_coeff=false)
	@assert j1 == 0 || j1 == 2
	@assert j2 == 0 || j2 == 2
	f0 = ell * (ell - 1) / ((2 * ell - 1) * (2 * ell + 1))
	f1 = - (2 * ell^2 + 2 * ell - 1) / ((2 * ell - 1) * (2 * ell + 3))
	f2 = (ell + 1) * (ell + 2) / ((2 * ell + 1) * (2 * ell + 3))
	if abs_coeff
		f0 = abs(f0)
		f1 = abs(f1)
		f2 = abs(f2)
	end
	if j1 == 0 && j2 == 0
		w = wldl[0,0]
	elseif j1 == 0 && j2 == 2
		w0 = wldl[0,-2]
		w1 = wldl[0, 0]
		w2 = wldl[0,+2]
		w = f0 * w0 + f1 * w1 + f2 * w2
	elseif j1 == 2 && j2 == 0
		w0 = wldl[-2,0]
		w1 = wldl[ 0,0]
		w2 = wldl[+2,0]
		w = f0 * w0 + f1 * w1 + f2 * w2
	else # j1 == 2 && j2 == 2
		w00 = wldl[-2,-2]
		w01 = wldl[-2, 0]
		w02 = wldl[-2,+2]
		w10 = wldl[ 0,-2]
		w11 = wldl[ 0, 0]
		w12 = wldl[ 0,+2]
		w20 = wldl[+2,-2]
		w21 = wldl[+2, 0]
		w22 = wldl[+2,+2]
		w = ( f0*f0*w00 + f0*f1*w01 + f0*f2*w02
		+ f1*f0*w10 + f1*f1*w11 + f1*f2*w12
		+ f2*f0*w20 + f2*f1*w21 + f2*f2*w22 )
	end
	return w
end


# convert wldl -> wljj
function wldl_to_wljj(infname, outfname; abs_coeff=false)
	wll = readdlm(infname)
	wldl = Dict()
	ell = Array{Int}(Int(sum(wll[:,end])))
	wjj = Array{Float64}(length(ell), 4)
	e0 = 2
	m2 = 3
	p2 = 4
	m4 = 5
	p4 = 6
	n = 1
	for m=1:size(wll,1)
		Bool(wll[m,end]) || continue  # if necessary Δℓ are not present, skip this
		ℓ = Int(wll[m,1])
		wldl[-2,-2] = wll[m-2,e0]
		wldl[ 0,-2] = wll[m  ,m2]
		wldl[ 2,-2] = wll[m+2,m4]
		wldl[-2, 0] = wll[m-2,p2]
		wldl[ 0, 0] = wll[m  ,e0]
		wldl[ 2, 0] = wll[m+2,m2]
		wldl[-2, 2] = wll[m-2,p4]
		wldl[ 0, 2] = wll[m  ,p2]
		wldl[ 2, 2] = wll[m+2,e0]
		ell[n] = ℓ
		wjj[n, 1] = wljj_dl(ℓ, 0, 0, wldl; abs_coeff=abs_coeff)
		wjj[n, 2] = wljj_dl(ℓ, 0, 2, wldl; abs_coeff=abs_coeff)
		wjj[n, 3] = wljj_dl(ℓ, 2, 0, wldl; abs_coeff=abs_coeff)
		wjj[n, 4] = wljj_dl(ℓ, 2, 2, wldl; abs_coeff=abs_coeff)
		n += 1
	end
	writedlm(outfname, [ell wjj])
	println("Created '$outfname'.")
end


# along ℓ
function calc_wlχR(χ, R, elllist=[2:100,112,125,150,200,300,400,500,600,700,800,900,1000,1200])
	calc_wldlχR(χ, R, elllist, 0, "data/wldl_chi$(χ)_R$(R).tsv")
	wldl_to_wljj("data/wldl_chi$(χ)_R$(R).tsv", "data/wljj_chi$(χ)_R$(R).tsv")
end


# along χ
function calc_wlRℓ(R, ℓ, χrange=logspace(0, 5, 800))
	calc_wldlRℓ(R, ℓ-2, 0, χrange; fname="data/wldl_R$(R)_ell$(ℓ-2).tsv")
	calc_wldlRℓ(R, ℓ,   0, χrange; fname="data/wldl_R$(R)_ell$(ℓ).tsv")
	calc_wldlRℓ(R, ℓ+2, 0, χrange; fname="data/wldl_R$(R)_ell$(ℓ+2).tsv")
	χ = readdlm("data/wldl_R$(R)_ell$(ℓ).tsv")[:,1]
	wll = Dict()
	wll[-2,-2] = readdlm("data/wldl_R$(R)_ell$(ℓ-2).tsv")[:,2]
	wll[-2, 0] = readdlm("data/wldl_R$(R)_ell$(ℓ-2).tsv")[:,4]
	wll[-2, 2] = readdlm("data/wldl_R$(R)_ell$(ℓ-2).tsv")[:,6]
	wll[ 0,-2] = readdlm("data/wldl_R$(R)_ell$(ℓ).tsv")[:,3]
	wll[ 0, 0] = readdlm("data/wldl_R$(R)_ell$(ℓ).tsv")[:,2]
	wll[ 0, 2] = readdlm("data/wldl_R$(R)_ell$(ℓ).tsv")[:,4]
	wll[+2,-2] = readdlm("data/wldl_R$(R)_ell$(ℓ+2).tsv")[:,5]
	wll[+2, 0] = readdlm("data/wldl_R$(R)_ell$(ℓ+2).tsv")[:,3]
	wll[+2, 2] = readdlm("data/wldl_R$(R)_ell$(ℓ+2).tsv")[:,2]
	w00 = wljj_dl(ℓ, 0, 0, wll)
	w02 = wljj_dl(ℓ, 0, 2, wll)
	w20 = wljj_dl(ℓ, 2, 0, wll)
	w22 = wljj_dl(ℓ, 2, 2, wll)
	writedlm("data/wljj_R$(R)_ell$(ℓ).tsv", [χ w00 w02 w20 w22])
end


# more tests, because some parameter combinations take a long time!
function calc_single()
	pk = PkSpectrum()
	v = wllrr(42, 42, 7158.442208311312, 6542.5979874801815, pk)  # fast

	v = wllrr(42, 42, 7158.442208311312, 6442.5979874801815, pk)  # slow
	println("v=$v, -1.7195370954759132e-9")

	v = wllrr(62, 58, 2303.0, 2533.3, pk)  # slow
	println("v=$v, -2.8238003384349807e-10")
end


end # module


WlLucas.calc_single()


WlLucas.calc_wlχR(2303.0, 1.1)
WlLucas.calc_wlχR(2303.0, 1.0, [2:1200])
WlLucas.calc_wlχR(2303.0, 0.9)
WlLucas.calc_wlχR(2303.0, 0.8)
WlLucas.calc_wlχR(2303.0, 0.7)
WlLucas.calc_wlχR(2303.0, 0.6)
WlLucas.calc_wlχR(2303.0, 0.5)
WlLucas.calc_wlχR(2303.0, 0.4)
WlLucas.calc_wlχR(2303.0, 0.3)
WlLucas.calc_wlχR(2303.0, 0.2)
WlLucas.calc_wlχR(2303.0, 0.1)

WlLucas.calc_wlRℓ(1.1, 42)
WlLucas.calc_wlRℓ(1.0, 42)
WlLucas.calc_wlRℓ(0.9, 42)
WlLucas.calc_wlRℓ(0.8, 42)
WlLucas.calc_wlRℓ(0.7, 42)
WlLucas.calc_wlRℓ(0.6, 42)
WlLucas.calc_wlRℓ(0.5, 42)
WlLucas.calc_wlRℓ(0.4, 42)
WlLucas.calc_wlRℓ(0.3, 42)
WlLucas.calc_wlRℓ(0.2, 42)
WlLucas.calc_wlRℓ(0.1, 42)

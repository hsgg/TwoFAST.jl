#!/usr/bin/env julia


module Leftovers


########### 2F1 Miller #####################

function Dellinv_mul_D!{T}(acinv::Complex{T}, dloc::T, fac::Complex{T}, z::T, D)
	const sqrt2 = sqrt(T(2))
	const oosqrt2 = 1 / sqrt2
	constA = sqrt2 * (acinv - dloc)
	constB = -oosqrt2 * acinv * z
	constC = (1 - dloc) * z
	Dnew11 = D[1,1] - constA * D[2,1]
	Dnew12 = D[1,2] - constA * D[2,2]
	Dnew21 = constB * D[1,1] + constC * D[2,1]
	Dnew22 = constB * D[1,2] + constC * D[2,2]
	D[1,1] = Dnew11 * fac
	D[1,2] = Dnew12 * fac
	D[2,1] = Dnew21 * fac
	D[2,2] = Dnew22 * fac
end


function calc_Dnellinv_back!{T}(ellbegin, ellend, R::T, a, dl, D)
	@assert ellbegin >= ellend

	z = R^2

	D[1,1] = 1
	D[1,2] = 0
	D[2,1] = 0
	D[2,2] = 1

	for ell=ellbegin-1:-1:ellend
		#println("Moving $(ell+1) -> $ell")
		c = ell + dl + 3/2
		acinv = a / c
		dloc = dl / c
		b = c + a - 1 - dl

		#Dl = [[1.0  -sqrt(2)*acinv]; [-acinv*z/sqrt(2)  z]]
		#D = Dl * D

		Dellinv_mul_D!(acinv, dloc, c / (R * b), z, D)
	end
end

function calc_Dnell_for!{T}(ellbegin, ellend, R::T, a, dl, D)
	@assert ellbegin <= ellend

	z = R^2

	D[1,1] = 1
	D[1,2] = 0
	D[2,1] = 0
	D[2,2] = 1

	for ell=ellbegin:ellend-1
		#println("Moving $(ell+1) -> $ell")
		c = ell + dl + T(3)/2
		acinv = a / c
		b = c + a - 1 - dl

		Dl = (R*b/c) * [
			[c/(c-a)*(c-dl)/(a+c-dl)     sqrt(T(2))*c/(c-a)*(a-dl)/(a+c-dl)/z];
			[a*c/(sqrt(T(2))*(c-a)*(a+c-dl))  c^2/((c-a)*(a+c-dl)*z)]]
		D[:] = Dl * D
	end
end


function calc_Pinf(oosqrt2=1/sqrt(2))
	Pinf = Array{Float64}(2,2)
	Pinf[1,1] = oosqrt2
	Pinf[1,2] = 1
	Pinf[2,1] = oosqrt2
	Pinf[2,2] = 0
	return Pinf
end

function calc_invPinf(sqrt2=sqrt(2))
	Pinf = Array{Float64}(2,2)
	Pinf[1,1] = 0
	Pinf[1,2] = sqrt2
	Pinf[2,1] = 1
	Pinf[2,2] = -1
	return Pinf
end


function adjust_matrix!(A)
	nx = 5
	ny = 5
	A[1,1] += nx*eps(real(A[1,1])) + ny*im*eps(imag(A[1,1]))
	A[2,2] += nx*eps(real(A[2,2])) + ny*im*eps(imag(A[2,2]))
	A[1,2] -= nx*eps(real(A[1,2])) + ny*im*eps(imag(A[1,2]))
	A[2,1] -= nx*eps(real(A[2,1])) + ny*im*eps(imag(A[2,1]))
end


const sqrt2 = sqrt(2)
const oosqrt2 = 1/sqrt(2)
const Pinf = calc_Pinf()
const invPinf = calc_invPinf()
function Pinftransform!(D)
	#D[:] = Pinf * D * inv(Pinf)
	oosqrt2D12 = oosqrt2 * D[1,2]
	D22 = D[2,2]
	D[1,2] = D[1,1] - oosqrt2D12 + sqrt2 * D[2,1] - D22
	D[2,1] = oosqrt2D12
	D[2,2] = D[1,1] - oosqrt2D12
	D[1,1] = oosqrt2D12 + D22
end


function calc_Anell(ellbegin, ellend, R, a, dl; direction=0, imax=100)
	# direction:
	#	-1 means higher down to lower ell
	#	0  means same direction as ellbegin -> ellend
	#	1  means lower up to higher ell
	reversed = false
	if (direction < 0 && ellbegin < ellend) || (direction > 0 && ellbegin > ellend)
		tmp = ellend
		ellend = ellbegin
		ellbegin = tmp
		reversed = true
	end

	D = Array{Complex{Float64}}(2,2)
	if ellbegin < ellend
		calc_Dnell_for!(ellbegin, ellend, R, a, dl, D)
	else
		calc_Dnellinv_back!(ellbegin, ellend, R, a, dl, D)
	end

	if reversed
		try
			D = inv(D)
		catch
			warn("Adjusting matrix, R=$R, a=$a, dl=$dl")
			warn("Adjusting matrix, ellbegin=$ellbegin, ellend=$ellend")
			warn("Adjusting matrix, direction=$direction")
			while imax > 0 && (det(D) == 0 || !all(isfinite(inv(D))))
				adjust_matrix!(D)
				imax -= 1
			end
			adjust_matrix!(D)
			D = inv(D)
		end
	end

	Pinftransform!(D)
	return D
end

function calc_Anell_fast(ellbegin, ellend, R, a, dl)
	D = Array{Complex{Float64}}(2,2)
	calc_Dnellinv_back!(ellbegin, ellend, R, a, dl, D)
	Pinftransform!(D)
	return D
end


function scale_seed!(fmatch, A, fseed)
	fmatch_seed = A * fseed
	lambda = fmatch ./ fmatch_seed
	fseed .*= lambda
end
function calc_fseed{T}(ellseed, ellmatch, fmatch, R::T, a, dl)
	#println("ellseed: $ellseed")
	A = calc_Anell_fast(ellseed, ellmatch, R, a, dl)
	#println("A: $A")

	if !all(isfinite.(A))
		# A not being finite means we hit the underflow gap.
		#println("A not finite (ellseed=$ellseed)")
		return [Complex{T}(Inf), Complex{T}(Inf)]
	end

	if det(A) == 0 || !all(isfinite.(inv(A)))
		# inv(A) having infinite elements does not necessarily mean we
		# hit the underflow gap! It just means the forward recursion is
		# unstable.
		fseed = [Complex{T}(1), Complex{T}(1)]
		#println("Using eigenvector (ellseed=$ellseed, fseed=$fseed)")
	else
		fseed = A \ fmatch
		#println("Inverting matrix (ellseed=$ellseed, fseed=$fseed)")
	end

	scale_seed!(fmatch, A, fseed)

	return fseed
end


function calc_ellseed{T}(R::T, ellmax::Int; acc=1e-10, pini=1e-14, prec=eps(T))
	@assert R < 1
	z = R^2
	if z^-ellmax < 1 / prec
		return 2ellmax + round(Int, (log(acc) - log(pini)) / log(z))
	else
		return ellmax + round(Int, (log(acc) - log(pini) + log(prec)) / log(z))
	end
end


function calc_fmax{T}(ellmax, ellmatch, fmatch, R::T, n, dl, alpha;
		fmax_tol::T=1e-10, imax=20000, growth=:exponential, elldiff=10)
	fmax = [Complex{T}(1), Complex{T}(1)]
	a = n / 2 + dl / 2
	if a == 0
                #println("a=0!")
		return fmax .* Mellell_pre(ellmax, ellmax+dl, R, n, alpha)
	elseif R == 1
                #println("R unity!")
		fell = calc_Mll_unity(ellmax, n, dl, alpha)
		return fell
	end
        #println("recurse!")
	rdiff = T(1)
	ellseed = ellmax
	#ellseed = calc_ellseed(R, ellmax) - elldiff
	while rdiff > fmax_tol && imax > 0
		imax -= 1
		ellseed += elldiff
		#println("ellseed: $ellseed")
		fseed = calc_fseed(ellseed, ellmatch, fmatch, R, a, dl)
		fmaxnew = calc_Anell_fast(ellseed, ellmax, R, a, dl) * fseed
		rdiff = norm(fmaxnew - fmax) / norm(fmaxnew)
		fmax[:] = fmaxnew
		if growth == :exponential
			elldiff = ceil(Int, 1.5*elldiff)
		elseif growth == :linear
			# keep elldiff = 10
		else
			error("unkown growth mode $growth")
		end
		#println("fseed:   $fseed")
		#println("fmaxnew: $fmaxnew")
		#println("rdiff: $rdiff")
		#exit(1)
		if !all(isfinite.(fseed))
			return fseed
		end
	end
	if imax == 0
		println("ellmax:   $ellmax")
		println("ellmatch: $ellmatch")
		println("fmatch: $fmatch")
		println("R: $R")
		println("n: $n")
		println("dl: $dl")
		println("fmax_tol: $fmax_tol")
		println("imax: $imax")
		error("Maximum iterations reached!")
	end
	return fmax
end


function calc_underflow_fmax{T}(ellmax, ellmatch, fmatch, R::T, n, dl, alpha)
	a = n / 2 + dl / 2
	lmin = 0
	lmax = ellmax
	lmid = div(lmin + lmax, 2)
	fmax = [Complex{T}(1), Complex{T}(1)]
	while lmid != lmin && lmid != lmax
		#println("==>$lmid:")
		#fseed = calc_fseed(lmid, ellmatch, fmatch, R, a, dl)
		fseed = calc_fmax(lmid, ellmatch, fmatch, R, n, dl, alpha;
			growth=:linear, elldiff=1)
		if !all(isfinite.(fseed)) || norm(fseed) < realmin(T)
			lmax = lmid
		else
			lmin = lmid
			fmax = fseed
		end
		#println("<==$lmid: $lmin<l<$lmax, $fseed")
		lmid = div(lmin + lmax, 2)
	end
	if maximum(abs.(fmax)) > 1e-100
		println("WARN at R=$R, n=$n, dl=$dl")
		println("WARN ellmax=$ellmax, ellmatch=$ellmatch")
		println("WARN fmatch=$fmatch")
		println("WARN lmin=$lmin, lmax=$lmax")
		println("WARN fmax = $fmax")
		println("WARN abs(fmax) = $(abs.(fmax))")
		warn("Assumption may be violated: Result is not close to zero")
		#fmax[1] = NaN
	end
	return fmax, lmin
end


function calc_2f1_RqmG_orig{T}(ell, R::T, dl; q=1.0, m::Int=500, G=log(1e4/1e-4), alpha=1e-4)
	t = 2 * T(pi) * m / G
	n = q - 1 - im * t
	a = n / 2 + dl / 2

	@assert R <= 1

	if R == 0 && ell+dl != 0
		return [Complex{T}(0), Complex{T}(0)], ell
	end

	B0 = Mellell_pre(0, 0+dl, R, n, alpha)
	f21 = calc_f0(R, n, dl)
	f0 = B0 * f21
	if !all(isfinite.(f0)) && n != 0 && R != 1
		warn("R: $R")
		warn("n: $n")
		warn("dl: $dl")
		warn("alpha: $alpha")
		warn("B0:  $B0")
		warn("f21: $f21")
		warn("f0:  $f0")
		error("Matchpoint could not be calculated")
	end
	fell = calc_fmax(ell, 0, f0, R, n, dl, alpha)
	#exit(1)
	#println("B0: $B0")
	#println("f0:   $f0")
	#println("fell: $fell")
	if !all(isfinite.(fell)) || norm(fell) < realmin(T)
		fell, ell = calc_underflow_fmax(ell, 0, f0, R, n, dl, alpha)
	end
	if !all(isfinite.(fell))
		Bl = Mellell_pre(ell, ell+dl, R, n, alpha)
		println("ERROR at ell=$ell, R=$R, q=$q, m=$m, G=$G, dl=$dl")
		println("B0: $B0")
		println("Bl: $Bl")
		println("f0:    $f0")
		println("farb0: $(B0.*Array{Complex{T}}(hyp2f1_bp1(0, dl, R, q, m, G)))")
		println("fell:    $fell")
		println("farbell: $(Bl.*Array{Complex{T}}(hyp2f1_bp1(ell, dl, R, q, m, G)))")
		error("Initial hypergeometric functions could not be calculated!")
	end
	return fell, ell
end

end

Leftovers.main()


# vim: set sw=4 et sts=4 :

module Miller

export miller


import Base.realmin
realmin(arr) = realmin(typeof(real(arr[1])))


# find_nseedmin():
#   This function finds an 'n' s.t. the approximatio n -> infty is valid.
#   We do this by testing that the eigenvalues are close to each other.
#   Specifically, we ensure that the eigenvalues are within an annulus in the
#   complex plane, where the radius of the annulus has a radius 'fracdist'
#   times the distance between the eigenvectors at infinity.
function find_nseedmin(laminf1, laminf2, BCfn; fracdist=1/2.1)
    rmax = fracdist * abs(laminf1 - laminf2)
    #println("laminf: $laminf1\t$laminf2")
    #println("rmax: $rmax")

    # test whether we are in the 'n->infty' limit:
    wearedone(n) = begin
        B, C, D, E = BCfn(n)
        lam1, lam2 = eigvals([B C; D E])
        #println("$n: $lam1\t$lam2")
        if abs(lam1 - laminf1) < rmax && abs(lam2 - laminf2) < rmax
            return true
        elseif abs(lam2 - laminf1) < rmax && abs(lam1 - laminf2) < rmax
            return true
        end
        return false
    end

    # increase 'n' until we reach a satisfactory point:
    n = 0
    while !wearedone(n)
        n = max(1, ceil(Int, 1.1n))
    end
    return n
end


# find_nseedmin():
#   This function finds a minimum nseed s.t. an error of order 1e16 dies out
#   when reaching 'nmax'.
function find_nseedmin2(nmax, BCfn)
    nseed = nmax
    r = 1.0
    while r > 1e-10
        nseed += 1
        B, C, D, E = BCfn(nseed)
        lam1, lam2 = eigvals([B C; D E])
        #println("$nseed($r): $lam1\t$lam2")
        al1 = abs(lam1)
        al2 = abs(lam2)
        if al1 > al2
            r *= al2 / al1
        else
            r *= al1 / al2
        end
    end
    return nseed
end


# estimate_prec_loss():
#   This function estimates the loss of precision by multiplying the ratio of
#   eigenvalues.
function estimate_prec_loss(nmax, BCfn)
    r = 1.0
    while nmax > 0
        B, C, D, E = BCfn(nmax)
        lam1, lam2 = eigvals([B C; D E])
        #println("$nseed($r): $lam1\t$lam2")
        al1 = abs(lam1)
        al2 = abs(lam2)
        if al1 > al2
            r *= al2 / al1
        else
            r *= al1 / al2
        end
        nmax -= 1
    end
    return r
end


# get_ndiffmin():
#   Calculate the number of iterations over which an error of order 1e16 dies
#   out.
function get_ndiffmin(laminf1, laminf2)
    al1 = abs(laminf1)
    al2 = abs(laminf2)
    if al1 > al2
        r = al2 / al1
    else
        r = al1 / al2
    end
    @assert r != 1
    #println("laminf: $laminf1\t$laminf2")
    prec = eps(typeof(real(laminf1)))
    dn = ceil(Int, log(prec) / log(r))
    return max(1, dn)
end


function calc_Amn_back(nbegin, nend, BCfn)
    @assert nbegin >= nend
    T = Complex{Float64}
    A = eye(T, 2)
    for n=nbegin-1:-1:nend
        B, C, D, E = BCfn(n)
        #An = [[B  C]; [D  E]]
        #A = An * A
        A11 = B * A[1,1] + C * A[2,1]
        A12 = B * A[1,2] + C * A[2,2]
        A21 = D * A[1,1] + E * A[2,1]
        A22 = D * A[1,2] + E * A[2,2]
        A[1,1] = A11
        A[1,2] = A12
        A[2,1] = A21
        A[2,2] = A22
    end
    return A
end


function scale_seed!(fmatch, A, fseed)
    fmatch_seed = A * fseed
    lambda = fmatch ./ fmatch_seed
    fseed .*= lambda
end


function calc_fseed{T}(nseed, BCfn, f0::T, fasymp::T)
    A = calc_Amn_back(nseed, 0, BCfn)
    #println("===> A: $A")

    if !all(isfinite.(A))
        # A not being finite means we hit the underflow gap.
        #println("===> A not finite (nseed=$nseed)")
        return T([Inf, Inf])
    end

    if det(A) == 0 || !all(isfinite.(inv(A)))
        # inv(A) having infinite elements does not necessarily mean we
        # hit the underflow gap! It just means the forward recursion is
        # unstable.
        #println("===> Returning asymptote")
        fseed = deepcopy(fasymp)
    else
        #println("===> estimating fmax")
        fseed = A \ f0
    end

    scale_seed!(f0, A, fseed)

    return fseed
end


function calc_fmax_fn{T}(nmax, BCfn, f0::T, fasymp::T, ndiffmin, nminseed;
                      fmax_tol=1e-10, imax=100, growth=:exponential, expfac=1.2)
    fmax = deepcopy(fasymp)
    rdiff = 1.0
    nseed = max(nmax, nminseed)
    ndiff = ndiffmin
    i = imax
    while rdiff > fmax_tol && i > 0
        imax -= 1
        nseed += ndiff
        fseed = calc_fseed(nseed, BCfn, f0, fasymp)
        fmaxnew = calc_Amn_back(nseed, nmax, BCfn) * fseed
        rdiff = norm(fmaxnew - fmax) / norm(fmaxnew)
        #println("nseed:   $nseed")
        #println("fseed:   $fseed")
        #println("fmax:    $fmax")
        #println("fmaxnew: $fmaxnew")
        #println("rdiff:   $rdiff")
        fmax[:] = fmaxnew
        if growth == :exponential
            ndiff = ceil(Int, expfac*(ndiff + 1 - ndiffmin)) + ndiffmin
        elseif growth != :linear
            error("unkown growth mode $growth")
        end
        if !all(isfinite.(fseed))
            return fseed
        end
    end
    if i == 0 && rdiff > fmax_tol
        warn("nmax:  $nmax")
        warn("nseed: $nseed")
        warn("f0: $f0")
        warn("fmax_tol: $fmax_tol")
        warn("rdiff: $rdiff")
        warn("imax,i: $imax,$i")
        warn("Maximum iterations reached!")
    end
    return fmax
end


function calc_underflow_fmax(nmax, calc_f)
    nmin = 0
    nmid = div(nmin + nmax, 2)
    # ensure we pass on the error when we find no solution:
    fmax = Array{Complex128}([NaN, NaN])
    while nmid != nmin && nmid != nmax
        #println("==>$nmid:")
        fnew = calc_f(nmid)
        if !all(isfinite.(fnew)) || norm(fnew) < realmin(fnew)
            nmax = nmid
        else
            nmin = nmid
            fmax = fnew
        end
        #println("<==$nmid: $nmin<n<$nmax, $fnew")
        nmid = div(nmin + nmax, 2)
    end
    if all(isfinite.(fmax)) && maximum(abs.(fmax)) > 1e-100
        warn("nmax=$nmax")
        warn("calc_f=$calc_f")
        warn("nmin=$nmin, nmax=$nmax")
        warn("fmax = $fmax")
        warn("abs(fmax) = $(abs.(fmax))")
        warn("Assumption may be violated: Result is not close to zero")
        return Array{Complex128}([NaN, NaN]), nmin
    end
    return fmax, nmin
end


function miller{T}(n, BCfn::Function, f0::T, fasymp::T, laminf1, laminf2;
                    calc_fmax=calc_fmax_fn, fmax_tol=1e-10)
    # try fast linear growth
    #println("trying fast linear growth:")
    #=
    precloss = estimate_prec_loss(n, BCfn)
    if precloss < 1e-3
        ndiffmin = get_ndiffmin(laminf1, laminf2)
    else
        ndiffmin = ceil(Int, n / 10)
    end
    =#
    ndiffmin = get_ndiffmin(laminf1, laminf2)
    #println("ndiffmin: $ndiffmin")
    fn1 = calc_fmax(n, BCfn, f0, fasymp, ndiffmin, 0, fmax_tol=fmax_tol,
                   imax=100, growth=:linear, expfac=1.5)
    all(isfinite.(fn1)) && return fn1, n

    # try linear growth at nmax
    #println("trying linear growth:")
    fn2 = calc_fmax(n, BCfn, f0, fasymp, 1, 0, fmax_tol=fmax_tol,
                   imax=100n+100ndiffmin, growth=:linear)
    all(isfinite.(fn2)) && return fn2, n

    # bisect
    #println("bisecting:")
    calc_f(n) = calc_fmax(n, BCfn, f0, fasymp, 1, 0; growth=:linear,
                          fmax_tol=fmax_tol, imax=100n+100ndiffmin)
    fn3, n3 = calc_underflow_fmax(n, calc_f)
    all(isfinite.(fn3)) && return fn3, n3

    warn("ERROR at n=$n")
    warn("n3: $n3")
    warn("ndiffmin: $ndiffmin")
    warn("BCfn: $BCfn")
    warn("f0:  $f0")
    warn("fn1: $fn1")
    warn("fn2: $fn2")
    warn("fn3: $fn3")
    warn("fasymp: $fasymp")
    warn("Initial value could not be calculated!")
    return f0, 0
end


end

# vim: set sw=4 et sts=4 :

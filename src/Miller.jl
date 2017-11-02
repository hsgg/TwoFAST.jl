module Miller

export calc_fn


import Base.realmin
realmin(arr::Array) = realmin(typeof(real(arr[1])))


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


function calc_fseed(nseed, BCfn, f0, fasymp)
    A = calc_Amn_back(nseed, 0, BCfn)

    if !all(isfinite.(A))
        # A not being finite means we hit the underflow gap.
        #println("A not finite (ellseed=$ellseed)")
        return typeof(f0)([Inf, Inf])
    end

    if det(A) == 0 || !all(isfinite.(inv(A)))
        # inv(A) having infinite elements does not necessarily mean we
        # hit the underflow gap! It just means the forward recursion is
        # unstable.
        fseed = deepcopy(fasymp)
    else
        fseed = A \ f0
    end

    scale_seed!(f0, A, fseed)

    return fseed
end


function calc_fmax_fn(nmax, BCfn, f0, fasymp; fmax_tol=1e-10, imax=20000,
                      growth=:exponential, ndiff=10)
    fmax = deepcopy(fasymp)
    rdiff = 1.0
    nseed = nmax
    while rdiff > fmax_tol && imax > 0
        imax -= 1
        nseed += ndiff
        fseed = calc_fseed(nseed, BCfn, f0, fasymp)
        fmaxnew = calc_Amn_back(nseed, nmax, BCfn) * fseed
        rdiff = norm(fmaxnew - fmax) / norm(fmaxnew)
        fmax[:] = fmaxnew
        if growth == :exponential
            ndiff = ceil(Int, 1.5*ndiff)
        elseif growth != :linear
            error("unkown growth mode $growth")
        end
        if !all(isfinite.(fseed))
            return fseed
        end
    end
    if imax == 0
        println("nmax:   $nmax")
        println("nmatch: $nmatch")
        println("fmatch: $fmatch")
        println("fmax_tol: $fmax_tol")
        println("imax: $imax")
        error("Maximum iterations reached!")
    end
    return fmax
end


function calc_underflow_fmax(nmax, calc_f, fasymp)
    nmin = 0
    nmid = div(nmin + nmax, 2)
    fmax = deepcopy(fasymp)
    while nmid != nmin && nmid != nmax
        #println("==>$nmid:")
        fseed = calc_f(nmid)
        if !all(isfinite.(fseed)) || norm(fseed) < realmin(fseed)
            nmax = nmid
        else
            nmin = nmid
            fmax = fseed
        end
        #println("<==$nmid: $nmin<n<$nmax, $fseed")
        nmid = div(nmin + nmax, 2)
    end
    if maximum(abs.(fmax)) > 1e-100
        println("WARN nmax=$nmax")
        println("WARN calc_f=$calc_f")
        println("WARN fasymp=$fasymp")
        println("WARN nmin=$nmin, nmax=$nmax")
        println("WARN fmax = $fmax")
        println("WARN abs(fmax) = $(abs.(fmax))")
        warn("Assumption may be violated: Result is not close to zero")
    end
    return fmax, nmin
end


function calc_fn(n, BCfn::Function, f0, fasymp; calc_fmax=calc_fmax_fn)
    fn = calc_fmax(n, BCfn, f0, fasymp)
    if !all(isfinite.(fn)) || norm(fn) < realmin(fn)
        calc_f(n) = calc_fmax(n, BCfn, f0, fasymp; growth=:linear, ndiff=1)
        fn, n = calc_underflow_fmax(n, calc_f, fasymp)
    end
    if !all(isfinite.(fn))
        println("ERROR at n=$n")
        println("BCfn: $BCfn")
        println("f0:     $f0")
        println("fn:     $fn")
        println("fasymp: $fasymp")
        error("Initial value could not be calculated!")
    end
    return fn, n
end


end

# vim: set sw=4 et sts=4 :

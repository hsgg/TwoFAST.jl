using TwoFAST
using QuadOsc
using Test

@testset "2-fast and quadosc comparison" begin

    # is log10(k) a problem?
    pk(k) = ...
    pk(log10_k) = ...


    function quadosc_cc_valid(ℓ, χ)
        a, E = quadosc(x->pk(log10(x))*Bessels.sphericalbesselj(ℓ, χ*x)^2 * (x^2), kmin, Inf, n->besselj_zero(ℓ+0.5, n)/χ)
        b, Eb = quadosc(x->pk(log10(x))*Bessels.sphericalbesselj(ℓ, χ*x)^2 * (x^2), kmax, Inf, n->besselj_zero(ℓ+0.5, n)/χ)
        return a
    end


    N = 4096
    chi0 = 0.001
    kmax = 1e4 #200/13 
    kmin = 4e-5 #2.5/7000
    q = 1.1
    ell = [2] 
    RR = [1.0]

    # calculate M_ll at high ell, result gets saved to a file:
    f21cache = F21EllCache(maximum(ell), RR, N; q=q, kmin=kmin, kmax=kmax, χ0=chi0)
    write("out/F21EllCache", f21cache)

    # calculate all M_ll, result gets saved to a file:
    mlcache = MlCache(ell, "out/F21EllCache", "out/MlCache")
    write("out/MlCache", mlcache)

    # calculate wljj:
    w00 = Array{Float64}(undef, N, length(RR))
    w02 = Array{Float64}(undef, N, length(RR))
    w22 = Array{Float64}(undef, N, length(RR))
    function outfunc(wjj, ell, rr, RR)
        w00[:,:] = wjj[1]
        #w02[:,:] = wjj[2]
        #w22[:,:] = wjj[3]
    end
    rr = calcwljj(pk, RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q, outfunc=outfunc, cachefile="out/MlCache/MlCache.bin");


end

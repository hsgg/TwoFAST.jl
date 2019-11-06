# TwoFAST

**Master Branch:**
[![Build Status](https://travis-ci.org/hsgg/TwoFAST.jl.svg?branch=master)](https://travis-ci.org/hsgg/TwoFAST.jl)
[![codecov.io](http://codecov.io/github/hsgg/TwoFAST.jl/coverage.svg?branch=master)](http://codecov.io/github/hsgg/TwoFAST.jl?branch=master)

**Next Branch:**
[![Build Status](https://travis-ci.org/hsgg/TwoFAST.jl.svg?branch=next)](https://travis-ci.org/hsgg/TwoFAST.jl)
[![codecov.io](http://codecov.io/github/hsgg/TwoFAST.jl/coverage.svg?branch=next)](http://codecov.io/github/hsgg/TwoFAST.jl?branch=next)

The 2-FAST (*2-point function from Fast and Accurate Spherical bessel
Transform*) algorithm is implemented here in the [Julia](https://julialang.org)
programming language.

The algorithm is documented in the paper [Fast and Accurate Computation of
Projected Two-point functions](https://arxiv.org/abs/1709.02401).

To install in Julia-1.0, press `]` to enter package mode, and then

```julia
   pkg> add TwoFAST
```


## Minimal example

Load the module:

```julia
    using TwoFAST
```

For both minimal examples we need a power spectrum. For example, we can use the
one in the `test/` subdirectory of this project:

```julia
    using Dierckx
    using DelimitedFiles
    path = homedir() * "/.julia/packages/TwoFAST/"
    path *= readdir(path)[1]
    path *= "/test/data/planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat"
    d = readdlm(path, comments=true)
    pk = Spline1D(d[:,1], d[:,2])
```

To calculate the real-space correlation function, use

```julia
    N = 1024    # number of points to use in the Fourier transform
    kmax = 1e3  # maximum k-value
    kmin = 1e-5 # minimum k-value
    r0 = 1e-3   # minimum r-value (should be ~1/kmax)

    print("ξ(r), ℓ=0, ν=0: ")
    r00, xi00 = xicalc(pk, 0, 0; N=N, kmin=kmin, kmax=kmax, r0=r0)

    print("ξ(r), ℓ=0, ν=-2:")
    r, xi0m2 = xicalc(pk, 0, -2; N=N, kmin=kmin, kmax=kmax, r0=r0)
```

To calculate the integrals over two spherical Bessel functions, we first
calculate the Fourier kernels at the highest needed ℓ. This is done with the
structure `F21EllCache`. Then, we generate the full *Mll*-cache for each ℓ.
This will automatically store the result in the file `out/MlCache/MlCache.bin`,
and all related info will be stored in the structure `MlCache`. Finally, to
actually calculate the *wljj*-terms we call the function `calcwljj()`. However,
to store the *wljj*-terms, we need to create the output arrays, and write a
function, `outfunc()`, that will store them in the arrays. The function
`outfunc()` will be called for each ℓ in the array `ell`. Here's an example:

```julia
    N = 4096
    chi0 = 1e-3
    kmin = 1e-5
    kmax = 1e3
    q = 1.1
    ell = [42]  # only ell=42 for this run
    RR = [0.6, 0.7, 0.8, 0.9, 1.0]

    # calculate M_ll at high ell, result gets saved to a file:
    f21cache = F21EllCache(maximum(ell), RR, N; q=q, kmin=kmin, kmax=kmax, χ0=chi0)
    write("out/F21EllCache", f21cache)

    # calculate all M_ll, result gets saved to a file:
    mlcache = MlCache(ell, "out/F21EllCache", "out/MlCache")
    write("out/MlCache", mlcache)

    # calculate wljj:
    w00 = Array{Float64}(undef, N, length(RR))
    w02 = Array{Float64}(undef, N, length(RR))
    function outfunc(wjj, ell, rr, RR)
        if ell == 42
            w00[:,:] = wjj[1]
            w02[:,:] = wjj[2]
        end
    end
    rr = calcwljj(pk, RR; ell=ell, kmin=kmin, kmax=kmax, N=N, r0=chi0, q=q, outfunc=outfunc, cachefile="out/MlCache/MlCache.bin")
```

# Tests

## Real-space projection

We compare with the `quadosc` algorithm.


## Spherical Bessel projection

There are several kinds of test:

    * Compare `TwoFAST` with previous versions of itself
                        
    * Compare `TwoFAST` with Nemo/Arb (useful to test 2F1 and recursion)
                        
    * Compare `TwoFAST` with Lucas~1995, both along ℓ and along χ, possibly
      along Δℓ

They have different levels of usefulness. Comparing TwoFAST with a previous
version of itself has the advantage of detecting algorithm changes. However,
the hardware also has a strong influence, so that this test is only useful if
developing on the machine that the tests were created on.

Comparing `TwoFAST` with Nemo/Arb is perfect for testing intermediate results,
such as the calculation of the Fourier kernels *Mll*, *2F1*, ... I will
definitely do this at some point. It is necessary for developing Miller's
algorithm in the Δℓ-direction. It is also useful for testing the recursion
relations.

Finally, comparing `TwoFAST` with other high-precision algorithms is useful for
testing the full chain from the power spectrum to the result.


### Along ℓ

The first test tests along ℓ. We repeat the test for several ratios R=χ'/χ and
several χ. Furthermore, we do the test for (j,j')=(0,0),(0,2),(2,0),(2,2), thus
we get some limited testing in the Δℓ direction.


### Along χ

Gotta do this! Mostly a matter of copying stuff...


### More derivatives

Also need to add tests for more derivatives, e.g. (j,j')=(4,4), etc.

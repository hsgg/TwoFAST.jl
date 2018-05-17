# Tests

There are several kinds of test:

    * Compare `TwoFAST` with previous versions of itself
                        
    * Compare `TwoFAST` with Nemo/Arb
                        
    * Compare `TwoFAST` with Quadosc, Lucas~1995, ...

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

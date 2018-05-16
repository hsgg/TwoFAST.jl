Julia Implementation
====================

The Julia implementation of the 2-FAST algorithm resides here. The code is
organized in two main modules, `Miller` and `TwoFAST`. Example code is in the
file `twofast_example.jl`. The module `Miller` is an implementation of, well,
Miller's algorithm, and the module `TwoFAST` contains the, uh, 2-FAST
algorithm. Next, we describe them in more detail, and also the module
`PerformanceStats`, which is a simple module for testing performance.


TwoFast
-------

The 2-FAST algorithm is split into three main tasks. First, the *Mll* are
calculated at *lmax* and *dlmax*. Second, the *Mll* (or rather the *Mljj*) are
calculated for all *l* and *dl*, and, finally, the *wljj* terms are calculated.

The implementation combines the first two in the sense that only one struct is
defined for them.

follows these three tasks by defining three structs that
define where the information is found if it was written to disk, whether to
calculate on-the-fly, etc. The default is to keep everything in memory.


Miller
------

The implementation of Miller's algorithm here is very specific to the 2-FAST
algorithm. Probably hard to adapt to anything else.


PerformanceStats
----------------

The module `PerformanceStats` provides the necessary machinery to use `@timed`
efficiently. That is, with this module you can use the `+=` operator on the
object returned by the `@timed` macro. That is, you can write the following
code:

    using PerformanceStats
    stats = @timed nothing   # initialize, if that's convenient
    stats += @timed res = somefunction()
    stats += @timed res = someotherfunction()
    timed_println("resources: ", stats)


Tests
-----

Finally, we describe some tests.

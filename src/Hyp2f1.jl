# Calculate the Gauss hypergeometric function with *sufficient* precision.

module Hyp2f1

export hyp2f1


############### Nemo ######################
#using Nemo
#CC = ComplexField(1024)
#toacb(x) = CC(real(x), imag(x))
#function hyp2f1(a::Number, b::Number, c::Number, z::Number)
#    #@show toacb(a) toacb(b) toacb(c) toacb(z)
#    res = Nemo.hyp2f1(toacb(a), toacb(b), toacb(c), toacb(z))
#    #@show res
#    return ComplexF64(Float64(real(res)), Float64(imag(res)))
#end


############### ArbNumerics ######################
using ArbNumerics
T = ArbComplex{1024}
function hyp2f1(a::Number, b::Number, c::Number, z::Number)
    res = ArbNumerics.hypergeometric_2F1(T(a), T(b), T(c), T(z))
    #return ComplexF64(Float64(real(res)), Float64(imag(res)))
    return Complex{Float64}(res)
end



end


# vim: set sw=4 et sts=4 :

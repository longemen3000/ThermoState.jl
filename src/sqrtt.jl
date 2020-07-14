
using BenchmarkTools, Random
Random.seed!(123)
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 1.0e-11 # may or may not matter
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1 # important

#= =
alternative implementation of square root.
#aproximation of the square root using the float representation as a trick
used in providing an initial value for calculations of square roots.
= =#

function aprox_sqrt1(z::Float64)
    val_int = reinterpret(Int64, z)
    val_int -= 1 << 52
    val_int >>= 1
    val_int += 1 << 61
    return reinterpret(Float64, val_int)
end
#= =
alternative implementation of square root.
#it uses a goldschmidt iteration scheme, with a sqrt aproximation as an initial value
five iterations on the goldschmidt part.
it takes advantage of the muladd instructions.
#the main problem is the slowness of the division instruction.
= =#
function goldschmidt_sqrt1(S::Float64)
    yn = 1.0 / aprox_sqrt1(S)
    xn = S * yn
    hn = yn / 2
    rn = 0.0
    for i = 1:5 # arbitrary number of loops
        rn = muladd(-xn, hn, 0.5)
        xn = muladd(xn, rn, xn)
        hn = muladd(hn, rn, hn)
    end
    return xn
end

#= =
Jeffrey Sarnoff:
Here is another sqrt developed from the continued fraction for sqrt 
and a fair amount of trial and error.
It takes 4/3 more time than goldschmidt_sqrt1; 
for the extra time (about 35ns vs 25ns on my machine)
one obtains the correctly rounded sqrt  for all rand(10_000) .* 1027 samples 
(goldschmidt_sqrt1 misses ~1/3 of them).

main cause for slowness? divisions?
= =#
function sqrtcf(x::Float64)
    y = aprox_sqrt1(x)  
    ry  = muladd(-y, y, x) / 2
    y1 = y + ry / y  
    ry  = muladd(-y1, y1, x) / 2
    y1 = y1 + ry / y1  
    ry  = muladd(-y1, y1, x) / 2
    a = y1 + ry / (2y1)
    b = y1 + ry / a
    return b
end
    

#= =one alternative implementation of inverse square root.
i cant make it work for now.
SIMPLE EFFECTIVE FAST INVERSE SQUARE ROOT ALGORITHM WITH TWO MAGIC CONSTANTS
Oleh Horyachyy, Leonid Moroz, Viktor Otenko
International Journal of Computing, 18(4) 2019, 461-470
ISSN 2312-5381
= =#
paper_inverse_square_root = "https://computingonline.net/computing/article/download/1616/883"


#= = 
john carmack inverse square root
the original code is for Float32. for Float64 the magic constant is other.

= =#
function fisr(x::Float64)
    x2 = 0.5 * x
    i = reinterpret(Int64, x)
    i = 6910469410427058089 - ( i >> 1 )
    y = reinterpret(Float64, i)
    y = y * ( 1.5 - ( x2 * y * y ) ) # first newton iteration
    # y = y * ( 1.5 - ( x2 * y * y ) ) #second newton iteration, you can repeat for more precision
    return y
end

#= = second alternative implementation of square root.
#it uses a goldschmidt iteration scheme, with the john carmack fast inverse square root as
a initial value. one newton iteration on the inverse square root and seven iterations on the 
goldschmidt part.
it takes advantage of the muladd instructions.
this last version is the fastest
= =#
function mysqrt2(S::Float64)
    x2 = 0.5 * S
    i = reinterpret(Int64, S)
    i = 6910469410427058089 - ( i >> 1 )
    yn = reinterpret(Float64, i)
    yn = yn * muladd(yn * yn, -x2, 1.5)
    # yn = yn * ( 1.5 - ( x2 * yn * yn )
    xn = S * yn
    hn = yn * 0.5
    rn = 0.0
    for _ in 1:7
        rn = muladd(-xn, hn, 0.5)
        xn = muladd(xn, rn, xn)
        hn = muladd(hn, rn, hn)
    end
    return xn
end
# ulp error
ulp_error(x, y) = abs(big(x) - big(y)) / big(eps(y))
#= =
David Sanders:
This will generate over the whole (huge) range of floats, including subnormals
excelent to test
= =#
N = 10^4
v = filter(!isnan, abs.(reinterpret.(Float64, rand(UInt64, N))))

#= = 
Test suite for correctness
#the base for comparison is sqrt(big(v))
#tested against sqrt(v)
= =#
function test_sqrt(f, v)
    println("test for " * string(f) * ":")
    my_v = f.(v)
    base_v = sqrt.(v)
    std_v = sqrt.(big.(v))
    my_error = abs.((my_v - std_v) ./ std_v)
    base_error = abs.((base_v - std_v) ./ std_v)
    @show maximum(my_error)
    @show maximum(base_error)
    println()
    @show sum(my_error) / N
    @show sum(base_error) / N
    println()
    @show maximum(ulp_error.(std_v, my_v))
    @show maximum(ulp_error.(std_v, base_v))
    println()
    return nothing
end

a = @benchmark sqrtcf.(t) setup = (t = v) # continued fraction
b = @benchmark mysqrt2.(t) setup = (t = v) # goldschmidt + inverse sqrt root
c = @benchmark sqrt.(t) setup = (t = v) # reference

test_sqrt(mysqrt2,v) # maximum error: 2.4e-16, mean error: 
test_sqrt(sqrtcf,v)


# a1 = 4602678819172646912 is the correct mask to eliminate exponents except 1
# a1 & number is the exponent

# -4602678819172646912 is the inverse
mexp(x) = reinterpret(Float64, reinterpret(Int64, x) & 4602678819172646912) # exponent
mexp2(x) = reinterpret(Float64, reinterpret(Int64, x) & ~(4607182418800017408)) # normalized number


#= =nearest iteration:
we have f = sqrt(x) 
and f0: 0 = f*f - x = muladd(f,f,-x)
it moves the result one eps if necesary
changing the last ifelse, you can have nearest, round up, or round down.    
= =#

function nearest_iteration(x)
    @assert isfinite(x) && !signbit(x)
    faithful_f = sqrt(x) # change here to a custom implementation of sqrt
    f0 = muladd(faithful_f, faithful_f, -x) # f0 for sqrt
    nearest_f = ifelse(iszero(f0), faithful_f, nextfloat(faithful_f, ifelse(signbit(f0), 1, -1)))
    return nearest_f
end

function roundup_iteration(x)
    @assert isfinite(x) && !signbit(x)
    faithful_f = sqrt(x) # change here to a custom implementation of sqrt
    f0 = muladd(faithful_f, faithful_f, -x) # f0 for sqrt
    roundup_f = ifelse(iszero(f0), faithful_f, nextfloat(faithful_f, ifelse(signbit(f0), 1, 0)))
    return roundup_f
end

function roundown_iteration(x)
    @assert isfinite(x) && !signbit(x)
    faithful_f = sqrt(x) # change here to a custom implementation of sqrt
    f0 = muladd(faithful_f, faithful_f, -x) # f0 for sqrt
    roundown_f = ifelse(iszero(f0), faithful_f, nextfloat(faithful_f, ifelse(signbit(f0), 0, 1)))
    return roundown_f
end

function LinearAlgebra.copy_oftype(A::AbstractArray{T,N}, ::Type{S}) where {T::Unitful.Quantity,N,} 
    return convert(AbstractArray{S,N}, A)
end

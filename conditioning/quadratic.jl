"""
Returns all solutions to ax^2 + bx + c.
Function will return NaN, if no real roots are found.
"""
function roots(a::Float64, b::Float64, c::Float64)
    @assert(c != zero(Float64), "c cannot be 0")
    r1 = zero(Float64)
    r2 = zero(Float64)
    n = num_roots(a, b, c)
    v = -b / (2 * a)
    println("vertex: ", v)
    if n == 0
        return NaN
    elseif n == 1
        return v
    else
        fv = a * v^2 + b * v + c
        x, y = initial_guess(a, b, c, v, fv)
        r1 = bisection(a, b, c, x, v)
        r2 = bisection(a, b, c, v, y)
    end

    return r1, r2
end

##
"""
Returns the number of real roots of the given quadratic
"""
function num_roots(a::Float64, b::Float64, c::Float64)
    d = b^2 - 4 * a * c
    if d < zero(Float64)
        return 0
    elseif d == zero(Float64)
        return 1
    else
        return 2
    end
end

##
"""
Computes two initial guesses for the bisection method.
"""
function initial_guess(a::Float64, b::Float64, c::Float64,
                       v::Float64, fv::Float64)
    # find left endpoint for left root
    x = v - 1
    fx = a * x^2 + b * x + c
    while sign(fv) == sign(fx)
        x -= 5
        fx = a * x^2 + b * x + c
    end

    # find right endpoint for right root
    y = v + 1
    fy = a * y^2 + b * y + c
    while sign(fv) == sign(fy)
        y -= 5
        fx = a * y^2 + b * y + c
    end

    return x, y
end

##
using Printf

"""
Computes the root located in the range of x and y using the Bisection method.
"""
function bisection(a::Float64, b::Float64, c::Float64, x::Float64, y::Float64)
    while x < (x + y) / 2 < y
        x, y = bisect_step(a, b, c, x, y)
    end

    return (x + y) / 2
end

function bisect_step(a::Float64, b::Float64, c::Float64, x::Float64, y::Float64)
    m = (x + y) / 2
    fx = a * x^2 + b * x + c
    fm = a * m^2 + b * m + c
    fy = a * y^2 + b * y + c
    if sign(fx) != sign(fm)
        return x, m
    else
        return m, y
    end
end

##
"""
Analytically computes the roots of a quadratic using the Citardauq Formula.
"""
function citardauq(a::Float64, b::Float64, c::Float64)
    num = 2 * c
    denom = -b + sqrt(b^2 - 4 * a * c)
    r1 = num / denom
    r2 = c / (a * r1)

    return r1, r2
end

##
#=
Here, we take an extremely poorly conditioned quadratic: x^2 + 200_000x - 10.
The true roots are very close to 0. We see that the bisection method takes an
extremely long time to converge to the true roots.
=#
a = 1.0
b = 200000.0
c = -10.0
b1, b2 = roots(a, b, c)
c1, c2 = citardauq(a, b, c)
@show b1, b2
@show c1, c2

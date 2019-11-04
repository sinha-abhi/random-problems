"""
Returns all solutions to ax^2 + bx + c.
Function will return NaN, if no real roots are found.
"""
function roots(a::Float64, b::Float64, c::Float64)
    @assert(c != zero(Float64), "c cannot be 0")
    n = num_roots(a, b, c)
    v = -b / (2 * a)
    if n == 0
        return NaN
    elseif n == 1
        return v
    else
        fv = a * v^2 + b * v + c
        x0, y0 = initial_guess(a, b, c, v, fv)
        println("x0: ", x0)
        println("y0: ", y0)
    end
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
    # find left root
    x = v - 5
    fx = a * x^2 + b * x + c
    while sign(fv) == sign(fx)
        x -= 5
        fx = a * x^2 + b * x + c
    end

    # find right root
    y = v + 5
    fy = a * y^2 + b * y + c
    while sign(fv) == sign(fy)
        y -= 5
        fx = a * y^2 + b * y + c
    end

    return x, y
end
## Testing
roots(1.0, 2.0, -10.0)

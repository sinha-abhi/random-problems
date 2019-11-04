"""
Returns all solutions to ax^2 + bx + c.
Function will return NaN, if no real roots are found.
"""
function roots(a::Float64, b::Float64, c::Float64)
    @assert(c != zero(Float64), "c cannot be 0")
    n = num_roots(a, b, c)
    if n == 0
        return NaN
    elseif n == 1
        return (-b / (2 * a))
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

function initial_guess(a::Float64, b::Float64, c::Float64)

end
## Testing
roots(-4.0, 12.0, -9.0)

"""
Exploration of the stability and accuracy of element-wise vector summation.
"""

"""
Simple method of element-wise sum. This algorith is indeed backwards stable.
"""
function naivesum(x::Vector{Float64})
    s = zero(Float64)
    for i = 1 : length(x)
        s += x[i]
    end
    return s
end

function kahan(x::Vector{Float64})
    s = zero(Float64)
    c = zero(Float64)
    for i = 1 : length(x)
            y = x[i] - c
            t = s + y
            c = (t - s) - y
            s = t
    end
    return s
end

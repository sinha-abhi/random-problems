module BuffonNeedle

"""
buffon-needle.jl
Abhi Sinha
October 07, 2019

Buffon's Needle implementation to estimate π (short needle).
In 1901, Mario Lazzarini used this method to come up with the approximation
π ≈ 355/113.

Example with distance between the lines of 1 and needle length of 1.
```
L = 1
len = 1
needle_prob = simulate(3408, L, len)
pi_est = (2 * len) / (L * needle_prob)
@show pi_est
```
"""

"""
We need a needle to drop. We can see fairly easily (from our assumptions) that
the y-coordinate of the needle doesn't actually matter.
   x       x-coordinate at the center of the needle
   θ       angle between the needle and the x-axis
   len     length of the needle
   x_left  "left" end of the needle
   x_right "right" end of the needle
"""
struct Needle
    x
    θ
    len
    x_left
    x_right
end

"""
Simluating a needle drop.
We can reduce the simulation to two vertical lines at x = 0 and x = L.
The x-coordinates of the ends of the needle can be computed as:
    x_left = x - len / 2 * sin(θ) and x_right = x + len / 2 * sin(θ).

function drop_needle
    len  length of needle
    L    distance between vertical lines
"""
function drop_needle(len, L)
    x = rand() * L
    θ = rand() * π
    x_left = x - len / 2 * sin(θ)
    x_right = x + len / 2 * sin(θ)
    return Needle(x, θ, len, x_left, x_right)
end

"""
Time to run our simulation.
A dropped needle intersected one of the vertical lines if
    a) x_left ≦ 0, or
    b) x_right ≧ L.

function simluation
    iter    number of iterations
    L       distance between vertical lines
    len     length of needle
"""
function simulate(iter, L, len)
    @assert(len <= L, "Needle length must be no more than the line distance")
    ninter = 0
    for i = 1 : iter
        needle = drop_needle(len, L)
        if needle.x_left <= 0 || needle.x_right >= L
            ninter += 1
        end
    end

    return ninter / iter
end

end # module

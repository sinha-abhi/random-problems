"""
Suppose you are standing at the center of an equilateral triangle with side
length 20 m. At each corner, there is a hungry raptor. The raptors are capable of
running at 15 m/s, except for one. He is wounded so he can only run at 10 m/s.
You run at 6 m/s. How long can you survive, and what direction should you run?

The original raptor problem: https://xkcd.com/135/.

Abhi Sinha, sinha45@purdue.edu
Jan 23, 2020
"""

using Base.Threads, LinearAlgebra, Printf

ate = 0.2                   # minimum distance before doom: 20 cm
r_dist = 20.0               # distance between the raptors
v = [6.0, 10.0, 15.0, 15.0] # velocites: [human, hurt raptor, raptor, raptor]
max_time = 3                # max simulation time
angle_change = 0.5          # time before angle change

function derivatives(angle, pos)
    dh = [cos(angle), sin(angle)] * v[1]
    dr = [(pos[1] - pos[i]) / norm(pos[1] - pos[i]) * v[i] for i in 2 : 4]

    return dh, dr
end

function single_step!(angle, dt, pos)
    dh, dr = derivatives(angle, pos)
    pos[1] += dh * dt
    pos[2] += dr[1] * dt
    pos[3] += dr[2] * dt
    pos[4] += dr[3] * dt
end

eaten(p) = norm(p[1]-p[2]) <= ate || norm(p[1]-p[3]) <= ate || norm(p[1]-p[4]) <= ate

#=
The layout of this simulation:
    R
        H   WR
    R
where R are healthy raptors, WR is the wounded raptor,
and H is the unlucky human. The location of the human is the origin.
=#
function simulate(angles, t_dir_change; niter::Int = 1000, verbose::Bool = false)
    time = 0.0
    dt = max_time / niter
    k = 1   # angle count
    angle = angles[k]
    pos = Array{Array}(undef, 4)

    # initial positions
    pos[1] = [0.0, 0.0]                         # human
    pos[2] = [1.0, 0.0] * r_dist                # hurt raptor
    pos[3] = [-0.5, sqrt(3.0) * 0.5] * r_dist   # raptor 1
    pos[4] = [pos[3][1], -pos[3][2]]            # raptor 2

    for i in 1 : niter
        single_step!(angle, dt, pos)
        time += dt
        if eaten(pos)
            verbose && @printf("Human caught in %f seconds\n", time)
            break
        end

        # direction change
        if time - floor(time) in t_dir_change
            k += 1
            try
                angle = angles[k]
            catch BoundsError
                # ignore
            end
        end
    end

    return time
end

# Simple test (ideal angle if direction change is not allowed)
println("Optimal straight line direction")
angles = [pi * 0.25]
t = simulate(angles, [0.5], verbose = true)

# What happens if we run toward the hurt raptor for half a second, and then
# go in the opposite direction?
println("Running toward the wounded raptor, and then the opposite direction")
angles = [0, pi]
simulate(angles, [0.5], verbose = true)

function grid_search_sampling(a, b, m)
    ranges = [range(a[i], stop = b[i], length = m[i]) for i in 1 : length(a)]
    collect.(collect(Iterators.product(ranges...)))
end

function raptor_grid_search(dims)
    a = Array{Int64}(undef, dims)
    b = Array{Int64}(undef, dims)
    m = Array{Int64}(undef, dims)

    best_time = Atomic{Float64}(-1.0)
    # do a grid search in batches of 5
    for i in 1 : 72 # 360 / 5
        batch_max = 5 * i
        fill!(a, batch_max - 5)
        fill!(b, batch_max - 1)
        fill!(m, 4)

        x = grid_search_sampling(a, b, m)
        @threads for s in x
            s = s .* (pi / 180) # convert to radians
            atomic_max!(best_time, simulate(s, [0.5], niter = 3000))
        end
    end

    return best_time
end

dims = convert(Int64, max_time / 0.5)
print(raptor_grid_search(dims))

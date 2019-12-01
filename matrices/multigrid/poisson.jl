"""
Solving Poisson's equation with Multigrid.
"""

include("multigrid_functions.jl")

using Plots
using SparseArrays

# Some basic benchmarking.
function time_poisson(n)
    nx = ny = n
    return @timed solve_poisson_direct(poisson_setup(nx, ny, (x, y) -> 1))
end

N = [2^i - 1 for i in 5 : 10]
results = map(n -> time_poisson(n), N)
solns = [results[i][1] for i in 1 : 6]
times = [results[i][2] for i in 1 : 6]

# TODO: plot solutions

plot(N, times, title = "Time to solve Poisson's Eq.", xaxis = :log, yaxis = :log)

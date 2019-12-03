"""
Solving Poisson's equation with Multigrid.
"""

include("multigrid_functions.jl")

using Plots
using SparseArrays

# Some basic benchmarking.
function time_direct_poisson(n)
    nx = ny = n
    return @timed solve_poisson_direct(poisson_setup(nx, ny, (x, y) -> 1))
end

N = [2^i - 1 for i in 5 : 10]
direct_results = map(n -> time__direct_poisson(n), N)
direct_solns = [direct_results[i][1] for i in 1 : 6]
direct_times = [direct_results[i][2] for i in 1 : 6]

# TODO: plot solutions

plot(N, direct_times, title = "Direct solve time", xaxis = :log, yaxis = :log)

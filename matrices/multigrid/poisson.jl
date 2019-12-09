"""
Solving Poisson's equation with Multigrid.
"""
include("multigrid_functions.jl")

using LinearAlgebra
using Plots
using SparseArrays

## Some basic benchmarking.
function time_direct_poisson(n)
    return @timed solve_poisson_direct(poisson_setup(n, n, (x, y) -> 1))
end

N = [2^i - 1 for i in 5 : 10]
nN = length(N)
direct_results = map(n -> time_direct_poisson(n), N)
direct_solns = [direct_results[i][1] for i = 1 : nN]
direct_times = [direct_results[i][2] for i = 1 : nN]

plot(N, direct_times, title = "Direct solve time", xaxis = :log, yaxis = :log)

## Error in interpolating a solution for Possion's eq.
F = poisson_setup(63, 63, (x, y) -> 1)
X = solve_poisson_direct(poisson_setup(31, 31, (x, y) -> 1))
Y = solve_poisson_direct(F) # true solution for 63x63
Y_est = interpolate(X) # interpolated solution for 63x63
E = abs.(Y - Y_est)
pyplot()
surface(E)
println("Norm of interpolation error: ", norm(E))

Y_est = apply_poisson_jacobi(Y_est, F)
E = abs.(Y - Y_est)
pyplot()
surface(E)
println("Norm of interpolation error after a Jacobi iteration: ", norm(E))

## Exploring simple_multigrid
niter = 25
_, E = simple_multigrid(31, 31, niter, false)
plot([1 : niter], E, title = "Error per iteration for n = 31")
_, E = simple_multigrid(63, 63, niter, false)
plot([1 : niter], E, title = "Error per iteration for n = 63")

## Solving via poisson_multigrid_residual(...)
function time_multigrid_residual(n, niter)
    return @timed poisson_multigrid_residual(n, n, niter)
end

multigrid_res = map(n -> time_multigrid_residual(n, 25), N)
mg_times = [multigrid_res[i][2] for i = 1 : nN]

plot(N, mg_times, title = "Multigrid solve time", xaxis = :log, yaxis = :log)

# Can we solve a very large system? - I can on my PC...
@time poisson_multigrid_residual(16383, 16383, 25)

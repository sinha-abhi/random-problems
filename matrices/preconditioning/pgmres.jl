"""
An exploration of Preconditioned GMRES.

Completed as part of David Gleich's Matrix Computations class found at
https://www.cs.purdue.edu/homes/dgleich/.
"""

## Load some data.
using DelimitedFiles
using SparseArrays

data = readdlm("poisson2D-data.csv", ',')
A = sparse(Int.(data[:, 1]), Int.(data[:, 2]), (data[:, 3]))
A = (A + A') ./ 2
b = vec(readdlm("poisson2D-rhs.csv"))

## Let's look at relative residuals of MINRES and GMRES without preconditioning.
using IterativeSolvers
using Plots

ε = 10^(-6) # convergence iff 2-norm of rel. residual is less than ε
maxit = 100
rest = 30
_, hist_min = minres(A, b, tol = ε, maxiter = maxit, log = true)
_, hist_gm = gmres(A, b, restart = rest, tol = ε, maxiter = maxit, log = true)
p = plot(1 : hist_gm.iters, hist_gm.data[:resnorm], label = ["GMRES"])
plot!(p, 1 : hist_min.iters, hist_min.data[:resnorm],
      label = ["MINRES"], linestyle = :dot, linewidth = 3)

# Now we work with some preconditioners.
include("preconditioners.jl")

M = LowerPrecond(ichol(A))
_, hist_pgm = gmres(A, b, restart = rest, tol = ε,
                    maxiter = maxit, Pl = M, log = true)

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
using Printf

ε = 10^(-6) # convergence iff 2-norm of rel. residual is less than ε
maxit = 100
rest = 30
_, hist_min = minres(A, b, tol = ε, maxiter = maxit, log = true)
_, hist_gm = gmres(A, b, restart = rest, tol = ε, maxiter = maxit, log = true)
p = plot(1 : hist_gm.iters, hist_gm.data[:resnorm], label = ["GMRES"])
plot!(p, 1 : hist_min.iters, hist_min.data[:resnorm],
      label = ["MINRES"], linestyle = :dot, linewidth = 3)

@printf("Iterations for MINRES: %d\n", hist_min.iters)
@printf("Iterations for GMRES: %d\n", hist_gm.iters)

# Now we work with some preconditioners.
include("preconditioners.jl")

# preconditioning MINRES
using Krylov
using LinearOperators

P = opInverse(LowerTriangular(ichol(A)))
Mp = P' * P
_, hist_pmin = Krylov.minres(A, b, M = Mp, itmax = maxit, rtol = ε)

# preconditioning GMRES
Mg = LowerPrecond(ichol(A))
_, hist_pgm = gmres(A, b, restart = rest, tol = ε,
                    maxiter = maxit, Pl = Mg, log = true)

@printf("Iterations for pMINRES: %d\n", length(hist_pmin.residuals))
@printf("Iterations for pGMRES: %d\n", hist_pgm.iters)

##
#=
A preconditioned system looks like: MAx = Mb. There are some interesting things
to see from the eigenvalues of A before and after preconditioning.
=#

A_m = Matrix(A) # It's fine to do this, since A isn't actually that big.
eig_A = eigvals(A_m)
plot(eig_A, zeros(length(eig_A)), seriestype = :scatter, title = "Eigvals of A")

# preconditioned matrix for MINRES
A_pm = Mp * A
eig_Apm = eigvals(Matrix(A_pm))
plot(eig_Apm, zeros(length(eig_Apm)), seriestype = :scatter, title = "Eigvals PMINRES")

# preconditioned matrix for GMRES
Mg = opInverse(LowerTriangular(ichol(A)))
A_pg = Mg * A
eig_Apg = [real(λ) for λ in eigvals(Matrix(A_pg))]
plot(eig_Apg, zeros(length(eig_Apg)), seriestype = :scatter, title = "Eigvals PGMRES")

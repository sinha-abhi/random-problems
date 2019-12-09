"""
Preconditioners formed by using Incomplete Cholesky and LU Factorizations.

Originally written by David Gleich (https://www.cs.purdue.edu/homes/dgleich/)
Some changes made by Abhi Sinha.
--
Example for forming preconditioners:

Lchol = LowerPrecond(ichol(Ahat)) # Cholesky Factorization
iLU = LUPrecond(ilu(Ahat)) # LU Factorization
"""

using LinearAlgebra

"""
Incomplete Cholesky factorization.
"""
function ichol(A::AbstractMatrix)
    L = copy(A)
    n = size(A, 1)

    for k = 1 : n
        L[k, k] = sqrt(L[k, k])
        for i = k + 1 : n
            if L[i, k] != 0
                L[i, k] /= L[k, k]
            end
        end

        for j = k + 1 : n
            for i = j : n
                if L[i, j] != 0
                    L[i, j] -= L[i, k] * L[j, k]
                end
            end
        end
    end
    tril!(L)

    return L
end

"""
Incomplete LU factorization.
"""
function ilu(B::AbstractMatrix)
    A = copy(B)
    n = size(A, 1)

    for i = 2 : n
        for k = 1 : i - 1
            if A[i, k] != 0
                A[i, k] /= A[k, k]
            end
            for j = k + 1 : n
                if A[i, j] != 0
                    A[i, j] -= A[i, k] * A[k, j]
                end
            end
        end
    end
    L = tril(A)
    U = triu(A)

    L -= Diagonal(L) + I

    return L,U
end

"""
Building the preconditioners.
"""
mutable struct LowerPrecond
    lower::LowerTriangular
end

mutable struct LUPrecond
    lower::LowerTriangular
    upper::UpperTriangular
end

"""
Some constructors.
"""
LowerPrecond(L::SparseMatrixCSC) = LowerPrecond(LowerTriangular(L))
LUPrecond(L::SparseMatrixCSC, U::SparseMatrixCSC) = LUPrecond(LowerTriangular(L), UpperTriangular(U))
LUPrecond(P::Tuple{SparseMatrixCSC{Float64, Int64}, SparseMatrixCSC{Float64, Int64}}) = LUPrecond(P[1], P[2])

"""
ldiv! implementations.
"""

import LinearAlgebra.ldiv!

ldiv!(L::LowerPrecond,x) = ldiv!(L.lower', ldiv!(L.lower,x))
ldiv!(y, L::LowerPrecond,x) = ldiv!(L, copy!(y,x))
ldiv!(P::LUPrecond,x) = ldiv!(P.upper, ldiv!(P.lower,x))
ldiv!(y, P::LUPrecond,x) = ldiv!(P, copy!(y,x))

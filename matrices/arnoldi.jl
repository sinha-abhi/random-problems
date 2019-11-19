using LinearAlgebra

"""
Searches the Krylov subspace of a symmetric matrix A of size nxn
the fact that λ_max = max{x^TAx}, where x is a unit eigenvector.
"""
function eigen_max(A::AbstractMatrix, b::AbstractVector)
    @assert(issymmetric(A), "A must be symmetric")
    @assert(size(A, 1) == size(b, 1), "vector b has incorrect dimensions")

    ε = 1e-12
    Q, H = arnoldi(A, b, ε)
    Λ = zeros(size(Q, 2))

    for i = 1 : size(Q, 2)
        q = Q[:, i]
        if norm(q - zeros(size(q))) < ε
            continue
        end
        α = q' * A * q
        println(α)
        Λ[i] = α
    end

    return maximum(Λ)
end

function arnoldi(A::AbstractMatrix, b::AbstractVector, ε)
    n = size(A, 1)
    Q = zeros(n, n + 1) # orthonormal basis for Krylov subspace
    H = zeros(n + 1, n) # basis for A on Q

    q = b / norm(b)
    Q[:, 1] = q
    for i = 1 : n
        q = A * q
        for k = 1 : i
            H[k, i] = Q[:, k]' * q
            q = q - H[k, i] * Q[:, k]
        end
        H[i + 1, i] = norm(q)
        if H[i + 1, i] > ε
            q = q / H[i + 1, i]
            Q[:, i + 1] = q
        else
            return Q, H
        end
    end

    return Q, H
end

## An example test case
A = [[1 2 3]; [2 3 4]; [3 4 5]]
b = [1; 1; 1]
Q, H = arnoldi(A, b, 1e-12)
λ = eigen_max(A, b)
α = eigmax(A)
println(λ)
print(α)
@show norm(α - λ)

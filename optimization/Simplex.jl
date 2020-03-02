"""
The Simplex Method

(Not the most efficient implementation of this algorithm, but it works.)

Adapted from David Gleich.
https://www.cs.purdue.edu/homes/dgleich/
"""

using LinearAlgebra

# Types
struct SimplexState
    c::Vector
    A::Matrix
    b::Vector
    bset::Vector{Int}
end

struct SimplexPoint
    x::Vector
    binds::Vector{Int}
    ninds::Vector{Int}
    lam::Vector
    sn::Vector
    B::Matrix
    N::Matrix
end

function SimplexPoint(T::Type)
    return SimplexPoint(zeros(T, 0), zeros(Int, 0), zeros(Int, 0), zeros(T, 0),
                        zeros(T, 0), zeros(T, 0), zeros(T, 0))
end

function SimplexPoint(T::Type, B::Matrix, N::Matrix)
    return SimplexPoint(zeros(T, 0), zeros(Int, 0), zeros(Int, 0),
                        zeros(T, 0), zeros(T, 0), B, N)
end

function simplex_point(state::SimplexState)
    m, n = size(state.A)
    @assert length(state.bset) == m "not enough indices to define a BFP"

    binds = state.bset
    ninds = setdiff(1 : size(A, 2), binds)
    B = state.A[:, binds]
    N = state.A[:, ninds]
    cb = state.c[binds]
    cn = state.c[ninds]
    c = state.c

    if rank(B) != m
        return (:Infeasible, SimplexPoint(eltype(c), B, N))
    end

    xb = B \ b
    x = zeros(eltype(xb), n)
    x[binds] = xb
    x[ninds] = zeros(eltype(xb), length(ninds))

    lam = B' \ cb
    sn = cn - N' * lam
    if any(xb .< 0)
        return (:Infeasible, SimplexPoint(x, binds, ninds, lam, sn, B, N))
    else
        if all(sn .>= 0)
            return (:Solution, SimplexPoint(x, binds, ninds, lam, sn, B, N))
        else
            return (:Feasible, SimplexPoint(x, binds, ninds, lam, sn, B, N))
        end
    end
end

function simplex_step!(state::SimplexState)
    # get the current point from the new basis
    stat, p::SimplexPoint = simplex_point(state)

    if stat == :Solution
        return stat, p
    elseif state == :Infeasible
        return :Breakdown, p
    else
        # Take the Dantzig index to add to basic
        qn = findmin(p.sn)[2]
        q = p.ninds[qn] # translate index
        # check that nothing went wrong
        @assert all(state.A[:, q] == p.N[:, qn])

        d = p.B \ state.A[:,q]

        # TODO: implement an anti-cycling method

        # Check for unbounded solutions
        if all(d .<= eps(eltype(d)))
            return :Degenerate, p
        end

        # determine the index to remove
        xq = p.x[p.binds] ./ d
        ninds = d .< eps(eltype(xq))
        xq[d .< eps(eltype(xq))] .= Inf
        pb = findmin(xq)[2]
        pind = p.binds[pb] # translate index

        @assert state.bset[pb] == pind
        state.bset[pb] = q

        return stat, p
    end
end

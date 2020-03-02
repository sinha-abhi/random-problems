"""
https://www.cs.purdue.edu/homes/dgleich/cs520-2020/julia/plotregion.jl
"""

module PlotRegion

using LazySets
using Combinatorics
using Plots
using LinearAlgebra

function basic_feasible_point(A::Matrix,b::Vector,set::Vector)
    m,n = size(A)
    @assert length(set) == m "need more indices to define a BFP"
    binds = set # basic variable indices
    ninds = setdiff(1:size(A,2),binds) # non-basic
    B = A[:,binds]
    N = A[:,ninds]
    #cb = c[binds]
    #cn = c[ninds]

    if rank(B) != m
        return (:Infeasible, 0)
    end

    xb = B\b
    x = zeros(eltype(xb),n)
    x[binds] = xb
    x[ninds] .= zero(eltype(xb))

    if any(xb .< 0)
        return (:Infeasible, x)
    else
        #lam = B'\cb
        #sn = cn - N'*lam
        return (:Feasible, x)
    end
end

"""
Plot the feasible polytope
Ax = b, x >= 0
for the first two components of x.
"""
function plotregion(A::Matrix,b::Vector)
    m,n = size(A)
    verts = Array{Array{Float64,1},1}(undef,0)
    for inds in combinations(1:n,m)
        bfp=basic_feasible_point(A,b,inds)
        if bfp[1] == :Feasible
            push!(verts,[bfp[2][1], bfp[2][2]])
        end
    end
    ch = convex_hull(verts)
    plot(VPolygon(ch),fillalpha=0.5, fillcolor="grey", label="")
end

export plotregion
end

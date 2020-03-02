include("Simplex.jl")
include("plotregion.jl")

using Plots
using Printf

function run_simplex(A, b, c, bset)
    PlotRegion.plotregion(A, b)

    state = SimplexState(c, A, b, bset)
    status, p = simplex_step!(state)
    iter = 1
    @printf("%d: Current p: (%f, %f)\n", iter, p.x[1], p.x[2])
    while status != :Solution && status != :Degenerate
        scatter!([p.x[1]], [p.x[2]], series_annotations = ["$(iter)"],
                 marker = (15, 0.2, :orange), label = "")
        status, p = simplex_step!(state)
        iter += 1
        @printf("%d: Current p: (%f, %f)\n", iter, p.x[1], p.x[2])
    end
    if status == :Degenerate
        println("WARNING: Degenerate solution")
    end
    @printf("Solution: (%f, %f)", p.x[1], p.x[2])
    scatter!([p.x[1]], [p.x[2]], series_annotations = ["$(iter)"],
             marker = (15, 0.2, :red), label = "")
end

#= Example 1
minimze         -5 x_1 - x_2
subject to      x_1 + x_2 <= 2
                2 x_2 + (1/2) x_2 <= 8
                x >= 0
=#
A1 = [1 1; 2 0.5]
b = [5; 8]
A = [A1 Matrix{Float64}(I, 2, 2)]
run_simplex(A, b, [-5, -1, 0, 0], [3, 4])

#= Example 2
minimize        -x_1 - 3 x_2
subject to      -2 x_1 + x_2 <= 2
                -x_1 + 2 x_2 <= 7
                x >= 0
=#
A1 = [-2 1; -1 2]
b = [2; 7]
A = [A1 Matrix{Float64}(I, 2, 2)]
run_simplex(A, b, [-1, -3, 0, 0], [3, 4])

#= Example 3
minimize        -(3/4) x_1 + 150 x_2 - (1/50) x_3 + 6 x_4
subject to      (1/4) x_1 - 60 x_2 - (1/25) x_3 + 9 x_4 <= 0
                (1/2) x_1 - 90 x_2 + (1/50) x_3 + 3 x_4 <= 0
                x_3 <= 1
                x >= 0
=#
A1 = [0.25 -60 -0.04 9;
      0.5  -90  0.02 3;
      0    0    1    0]
b = [0; 0; 1]
A = [A1 Matrix{Float64}(I, 3, 3)]
run_simplex(A, b, [-0.75, 150, -0.02, 6, 0, 0, 0], [5, 6, 7])

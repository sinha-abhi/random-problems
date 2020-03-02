"""
Booth function.
"""

using LinearAlgebra
using Optim
using Plots

booth(x) = (x[1] + 2 * x[2] - 7).^2 + (2 * x[1] + x[2] - 5).^2

# contour plot
pyplot()
x = -6. : 0.1 : 6.
y = -6. : 0.1 : 6.
X = repeat(x', length(y), 1)
Y = repeat(y, 1, length(y))
Z = map(booth, zip(X, Y))

surface(X, Y, Z)
contour(X, Y, Z)

optimize(booth, [0., 0.])
optimize(booth, [1000., 1000.])
optimize(booth, [3., 3.]) # 36 iterations

"""
Exponentiated Rosenbrock function.
"""

using LinearAlgebra
using Optim
using Plots

exprosen(x) = exp(100 * (x[2] - x[1].^2).^2 + (1 - x[1]).^2)
rlog(x) = log10(exprosen(x))

# semi-log contour plot
pyplot()
x = -2.5 : 0.1 : 2.5
y = -2.5 : 0.1 : 2.5
X = repeat(x', length(y), 1)
Y = repeat(y, 1, length(y))
Z = map(rlog, zip(X, Y))

surface(X, Y, Z)
contour(X, Y, Z)

# optimizating using the Optim package
optimize(exprosen, [0.0, 0.0])
optimize(exprosen, [-1.0, 1.0])
optimize(exprosen, [2.0, 2.0])

# Nelder-Mead exceeds max number of iterations
# What if try passing the gradient and Hessian?
function g!(grad, x)
    r = exprosen(x)
    grad[1] = (400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2) * r
    grad[2] = (-200 * x[1]^2 + 200 * x[2]) * r
end

function h!(hessian, x)
    grad = zeros(2)
    g!(grad, x)
    r = exprosen(x)

    _tmp = 400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2
    hessian[1, 1] = (1200 * x[1]^2 - 400 * x[2] + 2) * r + _tmp * grad[1]
    hessian[1, 2] = -400 * x[1] * r + _tmp * grad[2]

    _tmp = -200 * x[1]^2 + 200 * x[2]
    hessian[2, 1] = -400 * x[1] * r + _tmp * grad[1]
    hessian[2, 2] = 200 * r + _tmp * grad[2]
end

optimize(exprosen, g!, h!, [0.0, 0.0])

# Nope, Exponential Rosenbrock function near 3 is NaN
optimize(exprosen, g!, h!, [2.5, 2.5])
optimize(exprosen, g!, h!, [2.75, 2.75])
optimize(exprosen, g!, h!, [3.0, 3.0])

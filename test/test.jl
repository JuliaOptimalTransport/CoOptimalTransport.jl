using CoOptimalTransport
using Random
using Plots 
using LinearAlgebra
using Suppressor
using ImageTransformations

N = 10
M = 25

X1 = randn(N, M)
X2 = X1[randperm(size(X1, 1)), :][:, randperm(size(X1, 2))] # without replacement
X2 = hcat(X2, 2.5*randn(N, 5))

unif(n) = fill(1/n, n)
ε = 0.05
λ = 1.0

π_s, π_f = @suppress_err COOT(X1, X2, unif(size(X1, 1)), unif(size(X2, 1)), unif(size(X1, 2)), unif(size(X2, 2)), ε, ε; print_iter = 1);

plot(heatmap(π_s), heatmap(π_f))

π_s, π_f = @suppress_err UCOOT(X1, X2, unif(size(X1, 1)), unif(size(X2, 1)), unif(size(X1, 2)), unif(size(X2, 2)), ε, λ, λ; print_iter = 1);

plot(heatmap(π_s), heatmap(π_f))


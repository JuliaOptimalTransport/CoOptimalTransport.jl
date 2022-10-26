using Revise 
using CoOptimalTransport
using Random
using Plots 
using LinearAlgebra
using Suppressor
using ImageTransformations

N = 50
M = 25
X1 = randn(N, M)
X2 = X1 #[randperm(size(X1, 1)), :][:, randperm(size(X1, 2))] # without replacement
X2 = hcat(X2, 2randn(N, 10))
X2[rand(size(X2)...) .< 0.25] .= 0

unif(n) = fill(1/n, n)
ε = 0.05
λ = 1.0

π_s, π_f = @suppress_err COOT(X1, X2, unif(size(X1, 1)), unif(size(X2, 1)), unif(size(X1, 2)), unif(size(X2, 2)), ε, ε; print_iter = 1);
plt1=plot(heatmap(π_s), heatmap(π_f); plot_title = "COOT")
π_s, π_f = @suppress_err UCOOT(X1, X2, unif(size(X1, 1)), unif(size(X2, 1)), unif(size(X1, 2)), unif(size(X2, 2)), ε, λ, λ; print_iter = 1);
plt2=plot(heatmap(π_s), heatmap(π_f); plot_title = "UCOOT")

plot(plt1, plt2; layout = (2, 1))
savefig("COOT_vs_UCOOT.png")

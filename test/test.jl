using CoOptimalTransport
using Random
using Plots 

N = 10
M = 10

X1 = randn(N, M)
X2 = X1[randperm(size(X1, 1)), :][:, randperm(size(X1, 2))]

unif(n) = fill(1/n, n)

π_s, π_f = COOT(X1, X2, unif(size(X1, 1)), unif(size(X2, 1)), unif(size(X1, 2)), unif(size(X2, 2)), 0.05, 0.05);

heatmap(π_s)

sum(π_s)

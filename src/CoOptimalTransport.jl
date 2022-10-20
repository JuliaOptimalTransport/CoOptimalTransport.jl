module CoOptimalTransport

using OptimalTransport 
using LinearAlgebra
using LogExpFunctions

function COOT(X1, X2, w1, w2, v1, v2, εs, εf; iter = 100, print_iter = 25, tol = 1e-9)
    # initialize couplings
    π_s = w1 .* w2'
    π_f = v1 .* v2'
    π_s_new = similar(π_s)
    π_f_new = similar(π_f)
    T = float(Base.promote_eltype(X1, X2))
    Cs = Array{T}(undef, size(π_s)...)
    Cf = Array{T}(undef, size(π_f)...)
    for it = 1:iter
        # initialize cost matrices 
        Cs = -X1*π_f*X2'
        # OT for samples
        π_s_new .= sinkhorn(w1, w2, Cs, εs)
        Cf = -(X1')*π_s*X2
        π_f_new .= sinkhorn(v1, v2, Cf, εf)
        err = norm(π_s - π_s_new, 1) + norm(π_f - π_f_new, 1)
        if err < tol
            break
        end
        if it % print_iter == 0
            println("iteration $it, err = $err")
        end
        copy!(π_s, π_s_new)
        copy!(π_f, π_f_new)
    end
    return (π_s, π_f)
end

function UCOOT(X1, X2, w1, w2, v1, v2, ε, λ1, λ2; iter = 100, print_iter = 25, tol = 1e-9)
    # helper functions
    η = (x1, x2, p, q) -> (x1.^2 * p)/2 .+ (x2.^2 * q)'/2
    _sinkhorn_unbalanced(p, q, C, λ1, λ2, ε; kwargs...) = sinkhorn_unbalanced(p, q, C - ε*(1 .+ (log.(p) .+ log.(q)')), λ1, λ2, ε; kwargs...)
    H = (x, y) -> sum(xlogy.(x, x ./ y))
    # initialize couplings
    π_s = w1 .* w2'
    π_f = v1 .* v2'
    π_s_new = similar(π_s)
    π_f_new = similar(π_f)
    T = float(Base.promote_eltype(X1, X2))
    Cs = Array{T}(undef, size(π_s)...)
    Cf = Array{T}(undef, size(π_f)...)
    for it = 1:iter
        # π_f
        Cs = η(X1, X2, v1, v2) - X1*π_f*X2'
        Cs .+= ε*H(π_f, v1 .* v2') + λ1*H(vec(sum(π_f; dims = 2)), v1) + λ2*H(vec(sum(π_f; dims = 1)), v2)
        π_s_new .= _sinkhorn_unbalanced(w1, w2, Cs, λ1*sum(π_f), λ2*sum(π_f), ε*sum(π_f))
        # π_f
        Cf = η(X1', X2', w1, w2) - (X1')*π_s*X2
        Cf .+= ε*H(π_s, w1 .* w2') + λ1*H(vec(sum(π_s; dims = 2)), w1) + λ2*H(vec(sum(π_s; dims = 1)), w2)
        π_f_new .= _sinkhorn_unbalanced(v1, v2, Cf, λ1*sum(π_s), λ2*sum(π_s), ε*sum(π_s))
        err = norm(π_s - π_s_new, 1) + norm(π_f - π_f_new*sum(π_s_new)/sum(π_f_new), 1)
        if err < tol
            break
        end
        if it % print_iter == 0
            println("iteration $it, err = $err")
        end
        copy!(π_s, π_s_new)
        copy!(π_f, π_f_new*sum(π_s_new)/sum(π_f_new))
    end
    return (π_s, π_f)
end

export COOT
export UCOOT

end # module CoOptimalTransport

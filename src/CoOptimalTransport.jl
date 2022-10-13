module CoOptimalTransport

using OptimalTransport 

function COOT(X1, X2, w1, w2, v1, v2, εs, εf; iter = 100)
    # initialize couplings
    π_s = w1 .* w2'
    π_f = v1 .* v2'
    T = float(Base.promote_eltype(X1, X2))
    Cs = Array{T}(undef, size(π_s)...)
    Cf = Array{T}(undef, size(π_f)...)
    for it = 1:iter
        # initialize cost matrices 
        Cs = -X1*π_f*X2'
        # OT for samples
        π_s = sinkhorn(w1, w2, Cs, εs)
        Cf = -(X1')*π_s*X2
        π_f = sinkhorn(v1, v2, Cf, εf)
    end
    return (π_s, π_f)
end

export COOT

end # module CoOptimalTransport

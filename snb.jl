include("swf.jl")
function snb(G, s, l, r, eps)
    bds = 0
    n = G.n
    m = G.m
    for t=1:n
        if t!=s
            bds = bds + swf(G, s, t, l, r, eps)
        end
    end
    return bds
end

function snbp(G, s, l, r, phi, eps)
    n = G.n
    m = G.m
    theta = eps/2; tau = eps/2;
    p = phi * 1/sqrt(n) * sqrt(log(n))/tau
    bds = 0
    n = G.n
    m = G.m
    Vt = []
    for i in 1:n
        if rand() < p
            push!(Vt, i)  
        end
    end

    bds = 0
    for t in Vt
        if t!=s
            bds = bds + swf(G, s, t, l, r, theta)
        end
    end
    return bds/p
end
include("logw.jl")
include("graph.jl")
using Laplacians
using LinearAlgebra
using LightGraphs
using LightGraphs.SimpleGraphs
using LinearAlgebra
using SparseArrays
using Arpack
using Random
global_rng = Random.GLOBAL_RNG
rng = MersenneTwister(time_ns())


function push_bd(G, s, t, l)
    L = lapsp(G);
    A = adjsp(G);
    n = G.n
    m = G.m
    deg = zeros(n)
    for i=1:n
        deg[i] = L[i,i]
    end

    dm1 = spdiagm([1 / i for i in deg])
    pm = dm1 * A

    h = zeros(n,n);
    es = zeros(n); es[s] = 1; et = zeros(n); et[t] = 1;
    for i=0:l
        h += (es - et)
        es = es * pm
        et = et * pm
    end
    for i=1:n
        h[i] = h[i] / deg[i]
    end
    return norm(h)^2 - 1/n * sum(h)
end
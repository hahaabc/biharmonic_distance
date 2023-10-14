global_rng = Random.GLOBAL_RNG
rng = MersenneTwister(time_ns())
function random_walk(G, w, l, s)
    w[1] = s
    i = 1
    for ii = 0:l
        idx = rand(1:length(G.nbr[s]))
        s = G.nbr[s][idx]
        i = i + 1
        w[i] = s
    end
end

function swf(G, s, t, l, r, eps)
    
    L = lapsp(G);
    A = adjsp(G);
    n = G.n
    m = G.m
    deg = zeros(n)
    for i=1:n
        deg[i] = L[i,i]
    end
    S1 = zeros(n); S2 = zeros(n); T1 = zeros(n); T2 = zeros(n);
    Z = 0; sigma2 = 0; mk = 0; delta = 0.01;
    for k=1:r
        random_walk(G, S1, l, s)
        random_walk(G, S2, l, s)
        random_walk(G, T1, l, t)
        random_walk(G, T2, l, t)
        zk1 = 0; zk2 = 0;
        for x in S1
            for y in S1
                if x == y
                    if rand(rng) < 1/(deg[x])^2
                        zk1 += 1
                    end
                end
                if rand(rng) < 1/(deg[x] * deg[y])
                    zk2 += 1
                end
            end

            for y in T2
                if x == y
                    if rand(rng) < 1/(deg[x])^2
                        zk1 += 1
                    end
                end
                if rand(rng) < 1/(deg[x] * deg[y])
                    zk2 -= 1
                end
            end
        end

        for x in T1
            for y in T2
                if x == y
                    if rand(rng) < 1/(deg[x])^2
                        zk1 += 1
                    end
                end
                if rand(rng) < 1/(deg[x] * deg[y])
                    zk2 += 1
                end
            end

            for y in S2
                if x == y
                    if rand(rng) < 1/(deg[x])^2
                        zk1 += 1
                    end
                end
                if rand(rng) < 1/(deg[x] * deg[y])
                    zk2 -= 1
                end
            end
        end

        zk = zk1 - 1/n * zk2
        Z = Z + zk
        sigma2 = sigma2 + zk^2
        mk = Z/k
        sigma2k = sigma2^2/k - mk^2
        if sqrt(2*sigma2k * log(3/delta)/k) + 3 * psi * log(3/delta)/k < eps
            break
        end
    end
    return mk
end
global_rng = Random.GLOBAL_RNG
rng = MersenneTwister(time_ns())
function random_walk(G, s, i)
    for ii = 0:i-1
        idx = rand(1:length(G.nbr[s]))
        s = G.nbr[s][idx]
    end
    return s
end

function stw(G, s, t, l, r)
    for iiiii = 1:r
        f1 = 0
        f2 = 0
        f3 = 0
        f4 = 0
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        for i1 = 0:l
            for i2 = 0:l
                for rr = 1:r
                    s1 = random_walk(G, s, i1)
                    s2 = random_walk(G, s, i2)
                    if s1 == s2
                        if rand(rng) < 1/(deg[s1])^2
                            f1 += 1
                        end
                    end
                end
                for rr =1:r
                    t1 = random_walk(G, t, i1)
                    t2 = random_walk(G, t, i2)
                    if t1 == t2
                        if rand(rng) < 1/(deg[t1])^2
                            f2 += 1
                        end
                    end
                end

                for rr =1:r
                    t1 = random_walk(G, t, i1)
                    s2 = random_walk(G, s, i2)
                    if t1 == s2
                        if rand(rng) < 1/(deg[t1])^2
                            f3 += 1
                        end
                    end
                end    

                for rr =1:r
                    s1 = random_walk(G, s, i1)
                    t2 = random_walk(G, t, i2)
                    if s1 == t2
                        if rand(rng) < 1/(deg[s1])^2
                            f4 += 1
                        end
                    end
                end     

                for rr =1:r
                    s1 = random_walk(G, s, i1)
                    s2 = random_walk(G, s, i2)
                    t1 = random_walk(G, t, i1)
                    t2 = random_walk(G, t, i2)
                    if rand(rng) < 1/(deg[s1]*deg[s2])
                        n1 += 1
                    end
                    if rand(rng) < 1/(deg[t1]*deg[t2])
                        n2 += 1
                    end
                    if rand(rng) < 1/(deg[s1]*deg[t2])
                        n3 += 1
                    end
                    if rand(rng) < 1/(deg[t1]*deg[s2])
                        n4 += 1
                    end
                end
            end
        end
    end    
    println(s," ", t, " ", (f1+f2-f3-f4)/r - (n1+n2-n3-n4)/n/r)
    return (f1+f2-f3-f4)/r - (n1+n2-n3-n4)/n/r
end

    
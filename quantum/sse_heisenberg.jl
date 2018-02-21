#=  Implementation of SSE algorithm on Heisenberg model
    The code is referred to Anders Sandvik's tutorial
    author: kv
=#

function init_config(L)
    #Init spin configuration
    N = L * L
    randspin() = (Int8)[1,-1][rand(1:2)]
    s = [randspin() for _ in 1:N]
    return s
end

function make_lattice(L)
    # Make bond sites
    nn = L * L
    nb = 2 * nn
    bsites = zeros(Int32, 2, nb)
    for y1 in 0: L-1
        for x1 in 0: L-1
            s = x1+y1*L+1
            x2 = mod(x1+1, L)
            y2 = y1
            bsites[1, s] = s
            bsites[2, s] = x2+(y2)*L+1
            x2 = x1
            y2 = mod(y1+1, L)
            bsites[1, s+nn] = s
            bsites[2, s+nn] = x2+(y2)*L+1
        end
    end
    return bsites
end

function diagonal_update(spins::Array{Int8}, 
                            opstring::Dict,
                            bsites::Array{Int32},
                            L, beta, mm, nh)
    #= Diagonal update
        return:
            nh(Int32): order of series
            opstring(Array): Operator string
    =#
    # number of sites
    nn = L*L
    # number of bonds
    nb = 2*nn
    prob = 0.5*beta*nb
    invprob = 1.0/prob

    for i in 0:mm-1
        op = opstring[i]
        if (op == 0)
            # TODO: CHECK
            b = min(floor(Int32, rand()*nb) + 1, nb)
            #b = min(rand(1:nb), nb)
            if (spins[bsites[1, b]] != spins[bsites[2, b]])
                # float type L-n
                fLn = convert(Float64, mm-nh)
                if (prob >= fLn || prob >= rand() * fLn)
                    opstring[i] = 2*b
                    nh += 1
                end
            end
        elseif (mod(op, 2) == 0)
            fLn = convert(Float64, mm-nh+1)
            p = invprob*fLn
            if (p >= 1.0 || p >= rand())
                opstring[i] = 0
                nh -= 1
            end
        else
            b = div(op, 2)
            spins[bsites[1,b]]= -spins[bsites[1,b]]
            spins[bsites[2,b]]= -spins[bsites[2,b]]
        end
    end

    return nh, opstring
end

function loop_update(spins::Array{Int8},
                    opstring::Dict,
                    #opstring::Array{Int32},
                    bsites::Array{Int32},
                    L, beta, mm, nh)
    nn = L*L

    # DO WE NEED CHANGE THESE INTO DICT?, Yes
    firstspinop = -ones(Int32, mm)
    lastspinop = -ones(Int32, mm)
    vertexlist=Dict() # X(v)
    for i in 0:4*mm-1 vertexlist[i]=0 end
    #println(vertexlist)

    # Link Vertex List
    # indexing should be double checked
    for v0 in 0:4:4*mm-1
        v0_4 = div(v0, 4)
        op=opstring[v0_4]
        if (op != 0)
            b=div(op, 2)
            s1=bsites[1, b]
            s2=bsites[2, b]
            v1=lastspinop[s1]
            v2=lastspinop[s2]
            if (v1 != -1)
                vertexlist[v1]=v0
                vertexlist[v0]=v1
            else
                firstspinop[s1]=v0
            end
            if (v2 != -1)
                vertexlist[v2]=v0+1
                vertexlist[v0+1]=v2
            else
                firstspinop[s2]=v0+1
            end
            lastspinop[s1]=v0+2
            lastspinop[s2]=v0+3
        else
            vertexlist[v0: v0+3] = 0
        end
    end

    # Connect links across the boundary
    for s1 in 1:nn
        v1 = firstspinop[s1]
        if (v1 != -1)
            v2 = lastspinop[s1]
            vertexlist[v2] = v1
            vertexlist[v1] = v2
        end
    end

    for v0 in 0:2:4*mm-1
        if (vertexlist[v0] < 1) continue end
        v1 = v0
        if (rand() < 0.5)
            while true
                v1_4 = div(v1, 4)
                # bitwise exclusive or, for bit shift
                opstring[v1_4] = xor(opstring[v1_4],1)
                vertexlist[v1] = -1
                v2 = xor(v1, 1)
                v1 = vertexlist[v2]
                vertexlist[v2] = -1
                if (v1 == v0) break end
            end
        else
            while true
                # bitwise exclusive or, for bit shift
                vertexlist[v1] = 0
                v2 = xor(v1, 1)
                v1 = vertexlist[v2]
                vertexlist[v2] = 0
                if (v1 == v0) break end
            end
        end
    end

    # Metropolis
    for i in 1:nn
        if (firstspinop[i] != -1)
            if (vertexlist[firstspinop[i]] == -1)
                spins[i] = -spins[i]
            end
        else
            if (rand() < 0.5)
                spins[i] = -spins[i]
            end
        end
    end
end

function main()
    println("2D S=1/2 Heisenberg SSE program")
    # Parameters
    L = 4
    # number of spins?
    nn = L * L
    beta = 1.0
    mm = 20
    nh = 0
    isteps = 100
    msteps = 10000
    nbins = 10

    spins = init_config(L)
    bsites = make_lattice(L)

    #opstring = zeros(Int32, mm)
    # For clearness, use dictionary as lookup list
    opstring = Dict()
    #initialize operator string
    for i in 0:mm-1 opstring[i]=0 end

    for i in 1:isteps
        nh, opstring = diagonal_update(spins, opstring, bsites, L, beta, mm, nh)
        loop_update(spins, opstring, bsites, L, beta, mm, nh)

        # Cut-off opstring
        new_mm = nh + div(nh, 3)
        if new_mm > mm
            println("Grow the operator string to ", new_mm)
            for i in mm:new_mm-1
                opstring[i] = 0
            end
            mm = new_mm
        end
    end

    println("Equilibrum: ", mm)

    for j in nbins:
        E = 0.0
        for i in 1:msteps
            nh, opstring = diagonal_update(spins, opstring, bsites, L, beta, mm, nh)
            loop_update(spins, opstring, bsites, L, beta, mm, nh)
            E = E-nh
        end
        E = E / (beta*nn*msteps) +0.5
        println("Energy per site: ", E)
    end
end

main()
using PyPlot

function initialize(L::Int64)
    s = [Int8(rand(0:1) == 0? -1:1) for i in 1:L^2]
    return reshape(s, L, L)
end

function update(temp::Float64,s,L::Int64)
    N = L*L
    for step in 1:N
        i, j=rand(1:L),rand(1:L)
        s′=s[i,j]
        ΔE=2*J*s′*(s[(i+L-2)%L+1,j]+s[(i+L)%L+1,j]+s[i,(j+L-2)%L+1]+s[i,(j+L)%L+1])
        if (rand() < exp(-ΔE/temp))
            s[i,j] = -s′
        end
    end
end

function calM(s::Array{Int8})
    return sum(s)
end

function calE(s::Array{Int8}, L)
    E=0.0
    for i in 1:L
        for j in 1:L
            s′=s[i,j]
            E-=J*s′*(s[(i+L-2)%L+1,j]+s[(i+L)%L+1,j]+s[i,(j+L-2)%L+1]+s[i,(j+L)%L+1])
        end
    end
    return E/4.0
end

L = 8::Int64
N = L*L
J = 1
EQU=40000
MCS=20000
nBin=10
M_all=zeros(nBin)
E_all=zeros(nBin)
M_T=[]
E_T=[]
temps=linspace(3, 0., 31)
s=initialize(L)
@show(s)

tic()
for temp in temps
    # Thermalization
    for i in 1:EQU
        update(temp,s,L)    
    end

    for bin in 1:nBin
        M_bin=0.0
        E_bin=0.0
        for i in 1:MCS
            update(temp,s,L)
            M=calM(s)
            E=calE(s,L)
            M_bin+=abs(M)/N
            E_bin+=E/N
        end
        M_all[bin]=M_bin/Float64(MCS)
        E_all[bin]=E_bin/Float64(MCS)
    end
    M_ave=mean(M_all)
    E_ave=mean(E_all)
    
    println("$temp   $M_ave $E_ave")
    append!(M_T,M_ave)
    append!(E_T,E_ave)
end
plot(temps,M_T,temps,E_T)
toc()
function init_state(L::Int64)
    s = [Int8(rand(0:1) == 0? -1: 1) for _ in 1:L^2]
    return reshape(s, L, L)
end

function cal_mag(s::Array{Int8})
    return Float64(sum(s))
end

function cal_eng(s::Array{Int8}, L::Int64)
    E = 0.0
    for i in 1:L
        for j in 1:L
            s_ = s[i, j]
            E -= s_*(s[(i+L-2)%L+1,j]+s[(i+L)%L+1,j]+s[i,(j+L-2)%L+1]+s[i,(j+L)%L+1])
        end
    end
    return E/4.0
end

function update(s::Array{Int8}, L::Int64, T::Float64, E::Float64, M::Float64)
    """ Update, take current E & M
        then return E and M
    """
    E_ = E
    M_ = M
    N = L^2
    for site in 1:N
        i, j = rand(1:L), rand(1:L)
        s_ = s[i, j]
        dE = 2 * s_ * (s[(i-2+L)%L+1,j]
                        + s[(i+L)%L+1,j]
                        + s[i,(j-2+L)%L+1]
                        + s[i,(j+L)%L+1])
        if (rand() < exp(-dE/T))
            s[i,j]= -s_
            # update only when accepted
            E_ += dE
            M_ -= 2s_
        end
    end
    return E_, M_
end
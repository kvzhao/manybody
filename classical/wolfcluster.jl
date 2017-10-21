
function neighbors(L::Int64, loc::CartesianIndex{2})
    """Get nearest neighbors' sites
        How about its performance?
    """
    i, j = loc[1], loc[2]
    return [CartesianIndex((i+L-2)%L + 1, j), CartesianIndex((i+L  )%L + 1, j),
            CartesianIndex(i, (j+L-2)%L + 1), CartesianIndex(i, (j+L  )%L + 1)]
end

function grow_cluster!(loc::CartesianIndex{2}, spins::Array{Int8, 2}, s0::Int8, L::Int64, T::Float64)
    """Connect the bonds by recursively grow the cluster
    """
    # Flip the selected spin
    spins[loc] = -spins[loc]
    # loop through neighboring site, get all neighboring sites
    for nbor in neighbors(L, loc)
        if spins[nbor] == s0 && rand() > exp(-2.0/T)
                grow_cluster!(nbor, spins, s0, L, T)
        end
    end
end

function ClusterUpdate!(spins::Array{Int8, 2}, iters, L::Int64, T::Float64)
    for i in 1:iters
        # Random start
        start = CartesianIndex((rand(1:L), rand(1:L)))
        s0 = spins[start]
        # not connected bonds
        grow_cluster!(start, spins, s0, L, T)
        # Calculate Energy change
    end
end

function flip_cluster!(spins::Array{Int8, 2}, cluster::BitArray{2})
    L = size(spins, 1)
    for j in 1:L, i in 1:L
        if cluster[i,j] spins[i,j] *= -1 end
    end
end

function changed_spins(spins::Array{Int8, 2}, cluster::BitArray{2})
    """
    """
    sum = 0
    L = size(spins,1)
    for j in 1:L, i in 1:L
        if cluster[i, j] sum += spins[i,j] end
    end
    return sum
end

function neighbors(L::Int64, loc::CartesianIndex{2})
    """Get nearest neighbors' sites
        How about its performance?
    """
    i, j = loc[1], loc[2]
    return [CartesianIndex((i+L-2)%L + 1, j), CartesianIndex((i+L  )%L + 1, j),
            CartesianIndex(i, (j+L-2)%L + 1), CartesianIndex(i, (j+L  )%L + 1)]
end

function grow_cluster!(loc::CartesianIndex{2}, spins::Array{Int8, 2}, cluster::BitArray{2}, L::Int64, T::Float64)
    """Connect the bonds by recursively grow the cluster
    """
    # Connect the bond, also we can directly flip the spin
    cluster[loc] = true
    # loop through neighboring site, get all neighboring sites
    for nbor in neighbors(L, loc)
        if cluster[nbor] == false && spins[nbor] == spins[loc] && rand() > exp(-2/T)
                grow_cluster!(nbor, spins, cluster, L, T)
        end
    end
end

function ClusterUpdate!(spins::Array{Int8, 2}, L::Int64, T::Float64)
    # Random start
    start = CartesianIndex((rand(1:L), rand(1:L)))
    # not connected bonds
    cluster = falses(size(spins))
    grow_cluster!(start, spins, cluster, L, T)
    # Calculate Energy change
    dE = -2*connected_spins(spins, cluster)
    if dE < 0 || rand() < exp(-dE) flip_cluster!(spins, cluster) end
end

function flip_cluster!(spins::Array{Int8, 2}, cluster::BitArray{2})
    L = size(spins, 1)
    for j in 1:L, i in 1:L
        if cluster[i,j] spins[i,j] *= -1 end
    end
end

function connected_spins(spins::Array{Int8, 2}, cluster::BitArray{2})
    sum = 0
    L = size(spins,1)
    for j in 1:L, i in 1:L
        if cluster[i, j] sum += spins[i,j] end
    end
    return sum
end

println(filter(x -> x > 4, [1,2,3,4,5,6,7,8]))

# condition, data
g = filter(x -> x > 0.5, [rand() for i in 1:10])
println(g)
fname = "outfile.dat"

A = randn(100, 100)
println(ndims(A))
println(size(A))

f = open(fname, "w")

println(f, A)

close(f)
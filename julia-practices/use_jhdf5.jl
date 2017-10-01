#Pkg.add("HDF5")
using HDF5

A = randn(1000, 1000)
println(size(A))

outfile = "data.h5"

h5open(outfile, "w") do file
    write(file, "A", A)
end

c = h5open(outfile, "r") do infile
    read(infile, "A")
end

println(size(c))
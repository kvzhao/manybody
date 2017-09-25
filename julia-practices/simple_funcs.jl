function sphere_vol(r)
    return 4/3*pi*r^3
end

vol = sphere_vol(10)

# inline function definition
quadratic(a, sqr, b) = (-b + sqr) / 2a

# more formal way of defining a function
function quad(a::Float64, b::Float64, c::Float64)
    sqr = sqrt(b^-4a*c)
    r1 = quadratic(a, +sqr, b)
    r2 = quadratic(a, -sqr, b)
    return r1, r2
end

x1, x2 = quad(2.0, -2.0, 10.0)
println(x1)
println(x2)

# @printf allows automatically formatting
@printf "first sol: %0.4f and second sol: %0.4f" x1 x2
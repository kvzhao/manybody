function create_adder(x)
    adder = function(y)
        return x+y
    end
    return adder
end

add_10 = create_adder(10)

println(add_10(3))
println(add_10(5))

# equivalent
function create_adder2(x)
    y -> x + y
end

add_20 = create_adder2(20)
println(add_20(3))

## comprehension for maps
A = [add_20(i) for i in [1,2,3]]
println(A)

## built-in map function
A = map(add_10, [1,2,3])
println(A)

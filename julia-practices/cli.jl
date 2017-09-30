println(length(ARGS))
if length(ARGS) == 1
    print(ARGS[1])
elseif length(ARGS) >= 2
    println(map(x->string(x, x), ARGS))
end
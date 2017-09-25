function showtypetree(T, level=0)
    println("\t" ^ level, T)
    for t in subtypes(T)
        if t != Any
            showtypetree(t, level+1)
        end
    end
end

showtypetree(Number)

@time begin
    nheads = @parallel (+) for i=1:2000000
        int(randbool())
    end
end

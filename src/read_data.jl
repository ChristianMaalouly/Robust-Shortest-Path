file_size = ["20", "40", "60", "80", "100", "120", "140", "160", "180", "200", "250", "300", "350", "400", "450", "500",
    "550", "600", "650", "700", "750", "800", "850", "900", "950", "1000", "1100", "1200", "1300", "1400", "1500",
    "1600", "1700", "1800", "1900", "2000", "2100", "2200", "2300", "2400", "2500"]
file_type = ["BAY", "COL", "NY"]

dataFiles = [s * "_USA-road-d." * t * ".gr" for s in file_size, t in file_type]

function getNeighbours(n, A)
    δplus = [[] for i = 1:n]
    δminus = [[] for i = 1:n]
    for (i, j) in A
        push!(δplus[i], j)
        push!(δminus[j], i)
    end
    return δplus, δminus
end

function read_data(file)
    f = open("../data/processed/" * file)
    n = parse(Int, strip(split(readline(f), " = ")[2]))
    s = parse(Int, strip(split(readline(f), " = ")[2]))
    t = parse(Int, strip(split(readline(f), " = ")[2]))
    S = parse(Int, strip(split(readline(f), " = ")[2]))
    d1 = parse(Int, strip(split(readline(f), " = ")[2]))
    d2 = parse(Int, strip(split(readline(f), " = ")[2]))
    p = parse.(Int, split(split(readline(f), " = ")[2][2:end-1], ","))
    ph = parse.(Int, split(split(readline(f), " = ")[2][2:end-1], ","))

    line = readline(f)
    line = readline(f)
    A = []
    d = zeros(Int, n, n)
    D = zeros(n, n)
    while line != ""
        i, j, w, v = split(line[1:end-1], " ")
        i = parse(Int, i)
        j = parse(Int, j)
        w = parse(Int, w)
        v = parse(Float64, v)
        d[i, j] = w
        D[i, j] = v
        push!(A, (i, j))
        line = readline(f)
    end
    close(f)
    return n, s, t, S, d1, d2, p, ph, A, d, D
end
using ProgressBars
using Plots
using ColorSchemes

dims = 100
steps = 300
A = ones(Float64, dims, dims)
A2 = ones(Float64, dims, dims)
B = zeros(Float64, dims, dims)
m = trunc(Int, dims/2)
B[m-1:m+1, 10:dims-10] .= 1.0
B[10:dims-10, m-1:m+1] .= 1.0
B2 = zeros(Float64, dims, dims)
cache = (A = [zeros(Float64, dims, dims) for _ in 1:steps], B = [zeros(Float64, dims, dims) for _ in 1:steps])

global feed = 0.055
global kill = 0.062
global D_A = 1.0
global D_B = 0.5
global delta = 1.0

function diffuse(coords::Tuple, C)
    
    #weights = [[0.05, 0.2, 0.05], [0.2, -1.0, 0.2], [0.05, 0.2, 0.05]]

    return sum( # this is not a perfect approximation
        get(C, (coords[1] + 1, coords[2]), C[coords[1], coords[2]]) *0.2 +
        get(C, (coords[1] - 1, coords[2]), C[coords[1], coords[2]]) *0.2 +
        get(C, (coords[1], coords[2] + 1), C[coords[1], coords[2]]) *0.2 +
        get(C, (coords[1], coords[2] - 1), C[coords[1], coords[2]]) *0.2 +
        get(C, (coords[1] + 1, coords[2] + 1), C[coords[1], coords[2]]) *0.05 +
        get(C, (coords[1] - 1, coords[2] + 1), C[coords[1], coords[2]]) *0.05 +
        get(C, (coords[1] - 1, coords[2] - 1), C[coords[1], coords[2]]) *0.05 +
        get(C, (coords[1] + 1, coords[2] - 1), C[coords[1], coords[2]]) *0.05 +
        get(C, coords, 1.0) * (-1.0)
    )

end

for i in ProgressBar(1:steps)

    global A
    global B

    Threads.@threads for x in 1:dims
        kill = 0.05 + x/dims * 0.02
        for y in 1:dims
            feed = 0.04 + y/dims * 0.04

            ABB = A[x, y]*B[x, y]*B[x, y]
            A2[x, y] = A[x, y] + (D_A * diffuse((x, y), A) - ABB + feed * (1 - A[x, y])) * delta
            B2[x, y] = B[x, y] + (D_B * diffuse((x, y), B) + ABB - (kill + feed) * B[x, y]) * delta

        end # for
    end # for

    cache.A[i] = deepcopy(A2)
    cache.B[i] = deepcopy(B2)

    A = cache.A[i]
    B = cache.B[i]

end # for

anim = @animate for i in ProgressBar(1:trunc(Int, steps/10))
    heatmap(cache.A[i*10], clim=(0,1), c = :gist_rainbow)
end # for

gif(anim, "reaction.gif", fps = 120)
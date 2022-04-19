using ProgressBars
using Plots
using ColorSchemes

theme(:dark)

# Constants
dims = 200
steps = 1000
interval = 40

# Variables for storing the state.
A = ones(Float64, dims, dims)
A2 = ones(Float64, dims, dims)
B = zeros(Float64, dims, dims)
B2 = zeros(Float64, dims, dims)

# Seed B
for _ in 1:dims
    B[max(1, trunc(Int, rand() * dims)), max(1, trunc(Int, rand() * dims))] = 0.9
end # for

# Setup caches
cache = (A = [zeros(Float64, dims, dims) for _ in 1:trunc(Int, steps/interval)], 
         B = [zeros(Float64, dims, dims) for _ in 1:trunc(Int, steps/interval)])

# Global constants
global feed = 0.055 
global kill = 0.062
global D_A = 1.0
global D_B = 0.5
global delta = 1.0

"""
    diffuse(coords, C) -> Float64

Returns the amount of chemical C (where `C` is a 2D matrix)
at `coords` after some diffusion.
"""
function diffuse(coords::Tuple, C)
        
    #weights = [[0.05, 0.2, 0.05], [0.2, -1.0, 0.2], [0.05, 0.2, 0.05]]

    #== 
    In order to calculate the new value for this cell:
    Create a weighted average of all the cells around it,
    (orthogonally adjacent cells have a weight of 0.2 and 
    diagonally adjacten cells have a weight of 0.05).
    Then subtract the value in the current cell.
    If any of the adjacent cells don't exist, assume
    they have the value of the current cell.
    ==#
    return sum(
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

# Looping in order to break up the task.
for loop in 1:20

    # Iterate for 'steps' many steps...
    for i in ProgressBar(1:steps)

        # Disambiguate access to A and B we defined above
        global A
        global B

        # For each x coordinate (columns)...
        Threads.@threads for x in 1:dims

            #kill = 0.05 + x/dims * 0.02
            
            # For each y coordinate (now rows)...
            # This therefore goes through every cell.
            for y in 1:dims

                #feed = 0.04 + y/dims * 0.04

                # Calculate the amount of ABB which will become BBB.
                ABB = A[x, y]*B[x, y]*B[x, y]

                #==
                Calculate the amount of A that should be at these coordinates
                and store it in A2 (so that all the new values are calculated
                before the old ones are used.
                ==#
                A2[x, y] = A[x, y] + (D_A * diffuse((x, y), A) - ABB 
                           + feed * (1 - A[x, y])) * delta

                # Do the same for B.
                B2[x, y] = B[x, y] + (D_B * diffuse((x, y), B) + ABB 
                           - (kill + feed) * B[x, y]) * delta

            end # for
        end # for

        # Store the current frame in the cache.
        # The chop is used to not store only necessary frames.
        chop = trunc(Int, i/interval - 1/(interval + 1)) + 1

        #== 
        Since arrays are just pointers in julia, a copy
        is created so that when A2/B2 is modified later it
        does not also modify the cache.
        ==#
        cache.A[chop] = deepcopy(A2)
        cache.B[chop] = deepcopy(B2)

        # Update A and B with their new values so that 
        # the next iteration can run.
        A = cache.A[chop]
        B = cache.B[chop]

    end # for

    # For each relevant frame in the cache...
    anim = @animate for frame in ProgressBar(cache.A)

        # Create a heatmap showing the amount of A.
        heatmap(frame, clim=(0,1), c = :gist_rainbow)
        
    end # for

    # Great a gif of all the heatmaps.
    gif(anim, "reaction"*string(loop)*".gif", fps = 30)

end # for
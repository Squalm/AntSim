using Plots
using ProgressBars

steps = 200
dims = 100
ants = 200
wander = 0.5
speed = 1
global dirs = [-pi/4, 0, pi/4]

ph = zeros(Float64, dims, dims)
ant = Vector{Float64}[Float64[dims/2, dims/2, rand()*2*pi] for _ in 1:ants]
prev_ant = Vector{Float64}[Float64[dims/2, dims/2, rand()*2*pi] for _ in 1:ants]
holding_food = fill(false, ants)
food = zeros(Float64, dims, dims)
space = fill(true, dims, dims)
hive = [trunc(Int, dims/2), trunc(Int, dims/2)]

food_cache = [zeros(Float64, dims, dims) for _ in 1:steps]
ph_cache = [zeros(Float64, dims, dims) for _ in 1:steps]
ant_cache = [zeros(Float64, dims, dims) for _ in 1:steps]

# Patches of food
food[35:40, 35:40] .= 1.0
food[20:30, 80:90] .= 1.0
food[60:65, 60:65] .= 1.0

function walk(a::Vector{Float64}, ph, strength = 0.2)
    
    # Check in front, right and left
    ant[3] += dirs[findmax([ph[togrid(a[1] + cos(a[3] + dir) * speed), togrid(a[2] + sin(a[3] + dir) * speed)] 
        for dir in dirs])] * strength

    return a
end

function togrid(coord::Float64)
    return round(Int, coord)
end

function update(a::Vector{Float64})
    return [
        a[1] + cos(a[3]) * speed, 
        a[2] + sin(a[3]) * speed, 
        a[3]]
end

function grabfood(a::Int, holding_food)

    if !holding_food[a]
        # check around

    end # if

    return a
end

#==
for i in ProgressBar(1:steps)

    global prev_ant

    # reduce all pheremone levels
    if i % 50 == 0
        for x in 1:length(ph)
            ph[x] = max(ph[x] - 0.1, 0.0)
        end # for
    end # if

    _ant_alt = deepcopy(ant)

    for a in range(1,length(eachrow(ant)))

        # If next to food and not holding food, pick up food
        found_food = false
        if !holding_food[a]
            for d in go
                if get(food, (ant[a, 1] + d[1][1], ant[a, 2] + d[1][2]), 0.0) > 0.0
                    holding_food[a] = true
                    food[ant[a, 1] + d[1][1], ant[a, 2] + d[1][2]] -= 0.2
                    found_food = true
                    break
                end # if
            end # for
        end # if

        if holding_food[a]
            # If next to hive, drop food
            for d in go
                if [ant[a, 1] + d[1][1], ant[a, 2] + d[1][2]] == hive
                    holding_food[a] = false
                    break
                end # If
                # Lay pheremones around
                if get(space, (ant[a, 1] + d[1][1], ant[a, 2] + d[1][2]), false)
                    ph[ant[a,1] + d[1][1], ant[a,2] + d[1][2]] += 0.1
                end # if
            end # If

            # Lay pheremones on current squiare
            ph[ant[a, 1], ant[a, 2]] = min(ph[ant[a, 1], ant[a, 2]] + 0.1, 1.0)
            # If holding food, move to hive
            #dir_home = ((ant[a, 1] - hive[1])^2 + (ant[a, 2] - hive[2])^2)^(1/2)
            # deeply sad way of getting home
            if ant[a, 1] != hive[1] && 
                    !([ant[a, 1] + sign(hive[1] - ant[a, 1]), ant[a, 2]] in collect(eachrow(ant)))
                ant[a, 1] += sign(hive[1] - ant[a, 1])
            elseif ant[a, 2] != hive[2] && 
                    !([ant[a, 1], ant[a, 2] + sign(hive[2] - ant[a, 2])] in collect(eachrow(ant)))
                ant[a, 2] += sign(hive[2] - ant[a, 2])
            end # if
        else

            # look around at pheremone levels
            go = [d[1] => get(ph, (ant[a, 1] + d[1][1], ant[a, 2] + d[1][2]), -Inf) for d in go]
            # move in that direction
            sorted = all([d[2] <= 0.0 for d in go]) ? sort(go, by=x -> rand()) : sort(go, by=x -> x[2])
            for d in sorted
                global prev_ant
                if ant[a, 1] + d[1][1] < dims && ant[a, 1] + d[1][1] > 0 &&
                        ant[a, 2] + d[1][2] < dims && ant[a, 2] + d[1][2] > 0 && 
                        !([ant[a, 1] + d[1][1], ant[a, 2] + d[1][2]] in collect(eachrow(ant))) &&
                        collect(eachrow(prev_ant))[a] != [ant[a, 1] + d[1][1], ant[a, 2] + d[1][2]]
                    ant[a, 1] += d[1][1]
                    ant[a, 2] += d[1][2]
                    break
                end # if
            end # for
        end # if
    end # for

    prev_ant = deepcopy(_ant_alt)
    food_cache[i] = deepcopy(food)
    ph_cache[i] = deepcopy(ph)
    for a in eachrow(ant)
        ant_cache[i][a[1], a[2]] = true
    end # for

end # for
==#

println("Rendering...")
anim = @animate for i in ProgressBar(1:steps)
    heatmap(ant_cache[i], clim=(0,1))
end
 
gif(anim, "ants.gif", fps = 20)

anim1 = @animate for i in ProgressBar(1:steps)
    heatmap(ph_cache[i], clim=(0,1))
end

gif(anim1, "ph.gif", fps = 20)

anim2 = @animate for i in ProgressBar(1:steps)
    heatmap(food_cache[i], clim=(0,1))
end # for

gif(anim2, "food.gif", fps= 100)
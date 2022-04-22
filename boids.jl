using Plots
using ProgressBars

theme(:dark)

"""
Returns an array of all the boids within a radius r about centre c.
"""
function nearby(c, boids, r)
    return [(b, sqrt((c[1] - b[1])^2 + (c[2] - b[2])^2 )) for b in boids if sqrt((b[1] % w - c[1])^2 + (b[2] - c[2])^2) <= r]
end

function nearby(c, boids)
    return [(b, sqrt((c[1] - b[1])^2 + (c[2] - b[2])^2 )) for b in boids]
end # function

function s(a, b)
    @assert (size(a) == size(b))
    return sqrt(sum([(a[d] - b[d])^2 for d in 1:length(a)]))
end # function

"""
    keepin(boid, w, h, [edgeforce]) -> boid::Vector{Float64}

Returns the boid, with force applied to keep it inside a box
with width `w` and height `h`.
"""
function keepin(boid, w, h, edgeforce = 1)

    # Define a margin within which boids will 
    # be pushed away from the edge.  
    margin = 20

    #== X values ==#
    # If the boid is within the margin of the edge on the left...
    if boid[1] < margin

        # Push it right using a linear function getting
        # greater as it moves left.
        boid[3] += edgeforce * (margin - boid[1])
    
    # If the boid is within the margin of the edge on the right...
    elseif boid[1] > w-margin

        # Push it left using the same function.
        boid[3] -= edgeforce * (boid[1] - (w - margin))
    
    end # if

    #== Y values ==#
    # The same but substituting left for bottom,
    # and right for top.
    if boid[2] < margin
        boid[4] += edgeforce * (margin - boid[2])
    elseif boid[2] > h-margin
        boid[4] -= edgeforce * (boid[2] - (h - margin))
    end # if

    return boid

end

"""
    limit(boid, speedlimit) -> boid::Vector{Float64}

Returns the boid where the overall speed is capped 
at the speedlimit (direction is preserved).
"""
function limit(boid, speedlimit)
    
    # Calculate the 2D speed of the boid.
    speed = sqrt(boid[3]^2 + boid[4]^2)

    # If the speed is greater than the speedlimit...
    if speed > speedlimit

        # Modify the velocities so that the direction
        # is the same but the speed is the spedlimit.
        boid[3] = boid[3] * speedlimit/speed
        boid[4] = boid[4] * speedlimit/speed

    end # if

    return boid
    
end

"""
    separation(boid, [boidpush]) -> boid::Vector{Float64}

Returns the current boid with velocity changed to be pushed
away from the closest boid to it (with a 1/x relationship).  
"""
function separation(boid::Vector{Float64}, boidpush = 5.0)

    # Get all the nearby boids along with how far away they are.
    n = nearby(boid[1:2], boids, vision)

    #==
	We do not want the current boid to try and avoid itself,
	so I am careful to remove it. 
	==#
    n = [x for x in n if s(x[1][1:2], boid[1:2]) > 0]

    # If the current boid can see any boids at all...
    if length(n) > 0

        # Extract the distances from the list of boids and distances.
        near = [x[1] for x in n]

        # Find the location of the nearest boid.
        closest = near[findmin([x[2] for x in n])[2]][1:2]

        # Find the distance from the current boid to the nearest boid.
        distance = s(closest, boid[1:2])

        #== 
        Change the velocity of the current boid by using the trip
        relationship between side lengths and the 1/distance function
        multiplied by some constant (boidpush). 
        ==#
        boid[3] += boidpush * (boid[1] - closest[1]) / distance^2
        boid[4] += boidpush * (boid[2] - closest[2]) / distance^2

    end # if

    return boid

end

"""
    cohesion(boid, [centerpull]) -> boid::Vector{Float64}

Returns the current boid with velocity changed to be pulled
towards the centre of mass of boids it can see.
"""
function cohesion(boid::Vector{Float64}, centerpull = 0.08)

    # Initialise an empty variable that can hold [x, y] coordinates.
    center = Float64[0.0, 0.0]

    # Get the locations of the boids the current boid can see.
    near = [x[1] for x in nearby(boid[1:2], boids, vision)]

    #==
    Iterate through the list of nearby boids:
    For each boid, add the x and y coordinates to
    the x and y parts of center.
    ==#
    for b in near
        center[1] += b[1]
        center[2] += b[2]
    end # for

    # If there are any boids nearby...
    if length(near) > 0

        # Divide the center by the number of boids nearby, 
        # to give an average value.
        center ./= length(near)

        #==
        Take the difference between the centre and the current boid
        and multipy that by some constant. Then apply each
        component to to current boid's velocity.
        ==#
        boid[3] += (center[1] - boid[1]) * centerpull
        boid[4] += (center[2] - boid[2]) * centerpull
    
    end # if

    return boid

end

"""
    alignment(boid, [turnfactor]) -> boid::Vector{Float64}

Returns the current boid with velocity changed to move towards
the same speed and direction of the boids around it.
"""
function alignment(boid::Vector{Float64}, turnFactor = 0.1)

    # Initialise a variable that can hold the average values.
    avg = [0.0, 0.0]

    near = [x[1] for x in nearby(boid[1:2], boids, vision)]
    
    #==
    Iterate through the list of nearby boids:
    For each boid, add the V_x and V_y coordinates to
    the V_x and V_y parts of the average.
    ==#
    for b in near
        avg[1] += b[3]
        avg[2] += b[4]
    end # for

    # If there are any boids nearby...
    if length(near) > 0

        # Divide the sum by the number of boids nearby, 
        # to give an average value.
        avg ./= length(near)

        # Move towards the average speeds by some factor.
        boid[3] += avg[1] * turnFactor
        boid[4] += avg[2] * turnFactor

    end # if

    return boid

end

function update(boid::Vector{Float64}, delta)

	# Boid position is the current position (in components)
	# plus the velocity * some time interval.
    return [boid[1] + boid[3] * delta, boid[2] + boid[4] * delta, boid[3], boid[4]]

end

# Initialise constants
w = 200
h = 200
steps = 1000
vision = 15
delta = 0.2
n = 50
speedlimit = 10

# Create n boids in random locations including
# just outside of thet bounding box.
boids = Vector{Float64}[[rand()*3/2*w-w/4, rand()*3/2*h-h/4, (rand()-0.5) * 5, (rand()-0.5) * 5] for i in 1:n] # x,y,dx,dy 

anim = @animate for _ in ProgressBar(1:steps)

    global boids

    boids = cohesion.(boids)
    boids = separation.(boids)
    boids = alignment.(boids)
    boids = keepin.(boids, w, h)
    boids = limit.(boids, speedlimit)
        
    boids = update.(boids, delta)

    # Pull out the x and y values to be plotted.
    x = [boid[1] for boid in boids]
    y = [boid[2] for boid in boids]

    # Plot the positions of all the boids on a scatter graph.
    scatter(x, y, xlims=(0, w), ylims=(0, h))

end # for

gif(anim, "boids.gif", fps = 40)

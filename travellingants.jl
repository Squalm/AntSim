using ProgressBars
using Plots
using Combinatorics
using StatsBase

theme(:dark)

"""
	gen_points(n::Int, dims::Int) -> Vector{Vector{Float64}}
Returns `n` points in `dims` dimensions.
"""
function gen_points(n::Int, dims::Int)
	return [[rand() for _2 in 1:dims] for _ in 1:n]
end # function

"""
	s(a, b) -> Float64
Returns the euclidean distance between two n-dimensional points, `a` and `b`.
"""
function s(a, b)
	@assert (size(a) == size(b))
	return sum([(a[i] - b[i])^2 for i in 1:length(a)])^(0.5)
end # function

"""
	brute(points::Vector{Vector{Float64}}) -> (order = Vector{Int}, best_score = Float64)

Uses brute force to find the shortest route through `points`, returning to the starting point.  
Since it does not matter where we start, the first and last points are both 1.
"""
function brute(points::Vector{Vector{Float64}})

	# A variable which will store the distance of 
	# the shortest route the algorithm has found.
	best_score = Inf

	# A Vector which will store the best route.
	best = Int[]

	# For each possible order...
	# (assuming the order starts and ends at 1)
	for order in ProgressBar(permutations(2:length(points), length(points) - 1))

		# Calculate the distance of the order
		distance = route_s(points, order)

		# If the distance is shorter than the current best...
		if distance < best_score

			# Save the new distance and best order
			best_score = distance
			best = deepcopy(order)

		end # if

	end # for

	# Add the first and last point 1 (which was assumed)
	prepend!(best, 1)
	append!(best, 1)

	return (order = best, distance = best_score)

end # function

"""
	route(points, ph, n) -> Vector{Vector{Int64}}

Generate `n` routes through `points` where the pheremone intensities are `ph`.  
Uses weighted random sampling.
"""
function route(points::Vector{Vector{Float64}}, ph::Matrix{Float64}, n::Int)

	# Create a varibale which will hold the routes
	orders = Vector{Int64}[]

	# For each route that needs to be created...
	for _ in 1:n

		#== Get a route ==#
		# Create a variable which will temporarily store the route.
		order = Int[0 for _ in 2:length(points)]

		# Store the point where the ant currently is.
		current = 1

		# For each point the ant needs to travel through...
		for i in 1:length(points)-1

			# Store possible points the ant can go to
			# which is the list of points excluding the ones already visited.
			_Ps = [x for x in 2:length(points) if !(x in order)]
			
			# Store the Pheromone values at those point from the current point.
			_Phs = [ph[current, x-1] for x in 2:length(points) if !(x in order)]

			# Pick the next point randomly (weighted by the pheronome trails).
			p = sample(_Ps, pweights(_Phs), 1)[1]

			# Save that new point.
			order[i] = p

			# Move the ant to the new point.
			current = p

		end # for

		# Store the route in the list of routes.
		push!(orders, order)

	end # for

	return orders

end # function

"""
	route_s(points, orders) -> Vector{Float64}

Returns the distance of each route `orders` through `points`.
"""
function route_s(points::Vector{Vector{Float64}}, orders::Vector{Vector{Int}})

	return Float64[route_s(points, order) for order in orders]
	
end

"""
	route_s(points, order) -> Float64

Returns the distance of the route `order` through `points`
(in n-dimensions).
"""
function route_s(points::Vector{Vector{Float64}}, order::Vector{Int})

	# Create a variable to store the distance.
	distance = 0.0

	# For each point in the route (plus an end point)...
	for i in 1:length(order)+1

		#==
		Add the distance between the previous point
		(assumed starts at 1) and then next point to
		the distance variable.
		==#
		distance += s(points[get(order, i-1, 1)], points[get(order, i, 1)])

	end # for

	return distance

end # function

"""
	lay_ph(orders, distances, points, greed) -> Matrix{Float64}

Returns a matrix like ph which shows the new pheremones laid down in orders.
"""
function lay_ph(orders::Vector{Vector{Int}}, distances::Vector{Float64}, points::Vector{Vector{Float64}}, greed::Float64)

	# Create an empty matrix where the new pheromones can be stored.
	empty = fill(0.0, length(points), length(points) - 1)

	# For each order...
	for j in 1:length(orders)

		#== Update the pheremones ==#

		# Create a variable to store the previous point.
		prev = 1

		# For each point in the route...
		for i in orders[j]

			# Add pheromones to the path between the
			# previous point and the next point.
			empty[prev, i-1] += (1/distances[j]) ^greed

			# Save the current point as the previous point.
			prev = i
		
		end # for

	end # for

	return empty
	
end

"""
	orientation(p, q, r) -> Int

Returns `1` if the points go clockwise, `-1` if Anticlockwise
and `0` if they are colinear.
"""
function orientation(p, q, r)

	# This calculation gives a value that is determined
	# by the orientation (direction) of the points.
	dir = ((q[2] - p[2]) * (r[1] - q[1])) - ((q[1] - p[1]) * (r[2] - q[2]))

	# If the direction is greater than 0...
	if dir > 0

		# Then the points go clockwise.
		return 1

	# If the direction is less than 0...
	elseif dir < 0

		# Then the points go anticlockwise.
		return -1

	# Otherwise the direction is 0
	else

		# The points are colinear.
		return 0

	end # if

end # function

"""
	colinearOnSegment(p, q, r) -> bool

Given the points are colinear, checks if q lies on pr.
"""
function colinearOnSegment(p, q, r)

	# Knowing that q lies on the same line pr,
	# check if it is between them or not on the segment...
	if q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1]) && 
			q[2] <= max(p[2], r[2]) && q[2] >= min(p[2], r[2])
		return true
	end # if
	
	return false

end # function

"""
	ifintersect(p1, q1, p2, q2) -> bool

Checks if p1q1 and p2q2 intersect.
"""
function ifintersect(p1, q1, p2, q2)

	# Find the 4 relevant orientations
	o1 = orientation(p1, q1, p2)
	o2 = orientation(p1, q1, q2)
	o3 = orientation(p2, q2, p1)
	o4 = orientation(p2, q2, q1)

	# If the orientations of the two pairs are different
	# then the lines intersect.
	if (o1 != o2) && (o3 != o4)
		return true
	end # if

	# If any of the Special colinear cases...
	if (o1 == 0 && colinearOnSegment(p1, p2, q1)) ||
			(o2 == 0 && colinearOnSegment(p1, q2, q1)) ||
			(o3 == 0 && colinearOnSegment(p2, p1, q2)) ||
			(o4 == 0 && colinearOnSegment(p2, q1, q2))
		return true
	end # if

	# Otherwise the lines do not intersect.
	return false

end # function

"""
	decross(order, points, [limit]) -> order::Vector{Int}

Return an order were `limit` many changes were 
made to remove cross lines from the `order` through
`points`.
"""
function decross(order, points, limit = 100)

	# Create a variable we can test later to see 
	# if the algorithm made a change.
	success = false

	#==
	This function is recursive but it needs to 
	know to stop at some limit. This also allows
	us to stop the algorithm making too many changes.
	==#
	if limit > 0

		# Iterate through all the individual paths between points
		# in each route.
		for i in 0:length(order)
			for j in i+2:length(order)

					# Avoid a special case where the algorithm makes
					# the same change forever.
					if !(i == 0 && j == length(order))

						# Find the relevant 4 points.
						# The first line is p1q1.
						p1 = points[get(order, i, 1)]
						q1 = points[get(order, i+1, 1)]

						# The second line is p2q2.
						p2 = points[get(order, j, 1)]
						q2 = points[get(order, j+1, 1)]

						# If the two lines intersect...
						if ifintersect(p1, q1, p2, q2)
													
							# Remove the cross by reflecting the order
							# of points in the loop of the intersection. 
							order[i+1:j] = reverse(order[i+1:j])
							
							# Track that there was a success.
							success = true

							# break out of the loop.
							break
			
						end # if

					end # if

			end # for
			if success
				break
			end # if
		end # for

	end # if

	# If one pair of lines was uncrossed...
	if success

		#==
		Run the algorithm again on the new order.
		Since news crosses may have been introduced by
		the uncrossing.
		==#
		return decross(order, points, limit-1)

	# Otherwise...
	else

		# End the recursion and return the final order.
		return order

	end # if

end # function

"""
	ants(points, [iters, n, evaporation, greed, optimise]) -> (order, distance, ph)

Biomimicry with ant colonisation!  
Returns a route by sending `n` ants around in a loop.
"""
function ants(points::Vector{Vector{Float64}}, iters = 500, n = 20, evaporation = 0.1, greed = 2.0, optimise = true)

	# Create a matrix which will store the pheromone intensities.
	ph = fill(10.0/length(points), length(points), length(points) - 1)

	# Begin iterating...
	for _ in ProgressBar(1:iters)

		# Generate n routes through the points using the pheromone values.
		orders = route(points, ph, n)

		# If each route should be optimised slightly...
		if optimise

			# For each of the suggested orders...
			for order in orders

				# Run optimisation (decrosses most crossed lines)
				# on each order, limited to making k changes.
				order = decross(order, points, 3)

			end # for

		end # if

		# Calculate the distances of each route.
		distances = route_s(points, orders)

		# Place pheronmones on each route using the distances.
		ph += lay_ph(orders, distances, points, greed)

		# Allow some fraction of the pheromones to evaporate
		# The fraction is constant in all places.
		ph = (1-evaporation) .* ph

	end # for

	# Save the final pheronomes
	ph_out = deepcopy(ph)

	#== Find the preferred route ==#

	# An array to store the preferred route.
	best = Int[0 for _ in 2:length(points)]

	# A variable to store the current location of the
	# ant which follows the preferred route.
	current = 1

	# For each stop the ant must make...
	for i in 1:length(points)-1

		# Set the pheromones intensities of points the ant has already visited to -Infinity
		ph[current, :] = [x+1 in best ? -Inf : ph[current, x] for x in 1:length(ph[current, :])]
		
		# Find the point with the highest pheromone intensity.
		p = findmax(ph[current, :])[2] + 1

		# Save this point.
		best[i] = p

		# Move the ant to the new point.
		current = p

	end # for

	# Add 1s to the beginning and end of route (assumed).
	prepend!(best, 1)
	append!(best, 1)

	return (order = best, distance = route_s(points, best[2:end-1]), ph = ph_out)

end # function

"""
	draw(points, order, [f])

Plot the route `order` through `points`.
"""
function draw(points::Vector{Vector{Float64}}, order::Vector{Int}, f = "route.png")

	# Extract the x and y values for plotting.
	x = [points[i][1] for i in order]
	y = [points[i][2] for i in order]

	# Plot.
	plt = plot(x, y, xlims=(0, 1), ylims=(0, 1), marker=:cross, title=string(length(points)) * " points")

	# Save as a .png.
	png(plt, f)

end # function

function draw(ph::Matrix{Float64})

	plt = heatmap(ph)
	png(plt, "route_ph.png")

end # function

t = gen_points(18, 2)

a = ants(t, 20000, 50, 0.1, 2.0, true)
println("Ants solution: " * string(a.order) * " Distance: " * string(a.distance))
draw(t, a.order)

final_pass = decross(a.order[2:end-1], t, 1000)
append!(final_pass, 1)
prepend!(final_pass, 1)

draw(t, final_pass, "route_cross.png")
draw(a.ph)

#b = brute(t)
#println("Brute solution: " * string(b))
#draw(t, b.order, "route_brute.png")
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

	best_score = Inf
	best = Int[]

	for order in ProgressBar(permutations(2:length(points), length(points) - 1))

		distance = 0.0
		prev = 1
		for i in append!(order, 1)
			distance += s(points[prev], points[i])
			prev = i
		end # for

		if distance < best_score
			best_score = distance
			best = deepcopy(order)
		end # if

	end # for

	prepend!(best, 1)
	return (order = best, distance = best_score)

end # function

"""
Heck yeah Biomimicry with ant colonisation!  
"""
function ants(points::Vector{Vector{Float64}}, n = 10000)

	ph = Vector{Float64}[[0.1 for _2 in 2:length(points)] for _ in 1:length(points)]

	for _ in ProgressBar(1:n)

		# Get a route
		order = Int[0 for _ in 2:length(points)]
		current = 1
		for i in 1:length(points)-1
	
			_Ps = [x for x in 2:length(points) if !(x in order)]
			_Phs = [ph[current][x-1] for x in 2:length(points) if !(x in order)]

			p = sample(_Ps, pweights(_Phs), 1)[1]
			order[i] = p
			current = p
	
		end # for

		# Get the distance of the route
		distance = 0.0
		prev = 1
		for i in append!(order, 1)
			distance += s(points[prev], points[i])
			prev = i
		end # for

		# Update the pheremones
		prev = 1
		for i in order[1:end-1]
			ph[prev][i-1] += exp(-distance)
			prev = i
		end # for

	end # for

	# find the preferred route
	best = Int[0 for _ in 2:length(points)]

	distance = 0.0

	current = 1
	for i in 1:length(points)-1

		ph[current] = [x+1 in best ? -Inf : ph[current][x] for x in 1:length(ph[current])]
		
		p = findmax(ph[current])[2] + 1
		best[i] = p

		distance += s(points[current], points[p])

		current = p

	end # for
	distance += s(points[best[end]], points[1]) 

	prepend!(best, 1)
	append!(best, 1)
	return (order = best, distance = distance)

end # function

"""
	draw(points, order)
Plot the route `order` through `points`.
"""
function draw(points::Vector{Vector{Float64}}, order::Vector{Int})

	x = [points[i][1] for i in order]
	y = [points[i][2] for i in order]

	plt = plot(x, y, xlims=(0, 1), ylims=(0, 1), marker=:cross, title=string(length(points)) * " points")
	png(plt, "route.png")

end # function

t = gen_points(13, 2)
#draw(t, brute(t))

a = ants(t, 500000)
#b = brute(t)
println("Ants solution: " * string(a))
#println("Brute solution: " * string(b))
draw(t, a.order)
#draw(t, b.order)
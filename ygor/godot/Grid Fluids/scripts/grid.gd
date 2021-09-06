extends Node2D

var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(64, 64)
var tile_size 		= Vector2(grid_size.x/squares_qtd.x, grid_size.y/squares_qtd.y)
var show_vectors 	= false
var show_grid 		= false
var timer = 0.0

var rho = 1.0
var gravity = Vector2(0, -100) 
var sub_steps = 10 # random value, maybe be lowered for performance improvement
var MAX_VELOCITY = 1000

var Particle = preload("res://scenes/particle.tscn")
var Vector = preload("res://scenes/vector.tscn")

var grid_vectors = []
var particles = []

class VectorClass:
	var pressure: float
	var velocity: Vector2
	var pos: Vector2
	
	func check_vel(m_val):
		if velocity.x > m_val: velocity.x = m_val
		if velocity.x < -m_val: velocity.x = -m_val
	
		if velocity.y > m_val: velocity.y = m_val
		if velocity.y < -m_val: velocity.y = -m_val

func get_velocity(pos):
	return Vector2(pos.x, pos.y)

func get_pressure(pos):
	return (pos.x+pos.y)/100

func copy_vector(obj):
	var vec = VectorClass.new()
	vec.pressure = obj.pressure
	vec.velocity = obj.velocity
	vec.pos = obj.pos
	return vec

func copy_list(lista):
	var new = []
	for vec in lista:
		new.append(copy_vector(vec))
	return new

func copy_grid(grid):
	var new = []
	for line in grid:
		new.append(copy_list(line))
	return new

func _ready():
	for x in range(squares_qtd.y):
		grid_vectors.append([])
		for y in range(squares_qtd.x):
			# middle point
			var pos = Vector2((y+1) * tile_size.x + tile_size.x/2, (x+1) * tile_size.y + tile_size.y/2)
			var vector = Vector.instance()
			vector.pos = pos
			vector.velocity = get_velocity(pos)
			vector.pressure = get_pressure(pos) 
			grid_vectors[x].append(vector)
			$vector_visualizer.add_child(vector)
		
		var front = copy_vector(grid_vectors[x][0])
		var back = copy_vector(grid_vectors[x][-1])
		front.pos.x = tile_size.x/2
		back.pos.x += tile_size.x
		grid_vectors[x].push_front(front)
		grid_vectors[x].push_back(back)
	
	var front_list = []
	var back_list = []
	for x in grid_vectors[0]:
		var vec = copy_vector(x)
		vec.pos.y = tile_size.y/2
		front_list.append(vec)
	for x in grid_vectors[-1]:
		var vec = copy_vector(x)
		vec.pos.y += tile_size.y
		back_list.append(vec)
	grid_vectors.push_front(front_list)
	grid_vectors.push_back(back_list)

func bilinear_interpolation_vel(pos):
	var i = floor((pos.y-tile_size.y/2)/tile_size.y) + 1
	var j = floor((pos.x-tile_size.x/2)/tile_size.x) + 1
	
	var q11 = grid_vectors[i][j]
	var q12 = grid_vectors[i][j+1]
	var q21 = grid_vectors[i+1][j]
	var q22 = grid_vectors[i+1][j+1]
	
	var x1 = q11.pos.x
	var y1 = q11.pos.y
	
	var x2 = x1 + tile_size.x
	var y2 = y1 + tile_size.y
	
	return (q11.velocity * (x2 - pos.x) * (y2 - pos.y) +
		q21.velocity * (pos.x - x1) * (y2 - pos.y) +
		q12.velocity * (x2 - pos.x) * (pos.y - y1) +
		q22.velocity * (pos.x - x1) * (pos.y - y1)
		) / ((x2 - x1) * (y2 - y1))

func bilinear_interpolation_press(pos):	
	var i = floor((pos.y-tile_size.y/2)/tile_size.y) + 1
	var j = floor((pos.x-tile_size.x/2)/tile_size.x) + 1
	
	var q11 = grid_vectors[i][j]
	var q12 = grid_vectors[i][j+1]
	var q21 = grid_vectors[i+1][j]
	var q22 = grid_vectors[i+1][j+1]
	
	var x1 = q11.pos.x
	var y1 = q11.pos.y
	
	var x2 = x1 + tile_size.x
	var y2 = y1 + tile_size.y
	
	return (q11.pressure * (x2 - pos.x) * (y2 - pos.y) +
		q21.pressure * (pos.x - x1) * (y2 - pos.y) +
		q12.pressure * (x2 - pos.x) * (pos.y - y1) +
		q22.pressure * (pos.x - x1) * (pos.y - y1)
		) / ((x2 - x1) * (y2 - y1))

func add_particle():
	var pos = get_global_mouse_position()
	var new_particle = Particle.instance()
	new_particle.position = pos
	particles.append(new_particle)
	add_child(new_particle)
	print('particle added')

func update_boundary(w):
	# update vertical
	for x in range(squares_qtd.y+2):
		w[x][0].velocity = w[x][1].velocity
		w[x][0].pressure = w[x][1].pressure
		
		w[x][-1].velocity = w[x][-2].velocity
		w[x][-1].pressure = w[x][-2].pressure
		
	#update horizontal 
	for y in range(squares_qtd.x+2):
		w[0][y].velocity = w[1][y].velocity
		w[0][y].pressure = w[1][y].pressure
		
		w[-1][y].velocity = w[-2][y].velocity
		w[-1][y].pressure = w[-2][y].pressure
		
func external_forces():
	return gravity # + bouancy (+ mouse_reppelant_force)

func add_force(w, delta):
	for x in range(1, squares_qtd.y+1):
		for y in range(1, squares_qtd.x+1):
			w[x][y].velocity += delta * external_forces()

func advect(w, timestep):
	var s = timestep/sub_steps
	for x in range(1, squares_qtd.y+1):
		for y in range(1, squares_qtd.x+1):
			var pos = w[x][y].pos
			var vel = w[x][y].velocity
			for _s in range(sub_steps):
				pos -= s*vel
				if pos.x < 0: pos.x = 0
				if pos.y < 0: pos.y = 0 
				if pos.x > grid_size.x: pos.x = grid_size.x 
				if pos.y > grid_size.y: pos.y = grid_size.y
				vel = bilinear_interpolation_vel(pos) 
			w[x][y].pressure = bilinear_interpolation_press(pos)

func diffuse(_w):
	pass #Passthrough function because smoke has no viscosity

func poisson_solver(div, x0, tol):
	var max_it = 20

	for _i in range(max_it):
		var old = x0.duplicate(true)
		var accum = 0.0
		for x in range(1, squares_qtd.y+1):
			for y in range(1, squares_qtd.x+1):
				x0[x][y] = 0.25*(x0[x+1][y] + x0[x-1][y] + x0[x][y+1] + x0[x][y-1] + div[x][y].pressure)
				accum += abs(old[x][y]-x0[x][y])
		if accum < tol:
			break
	
	return x0

func divergent(w):
	#Refactor make func
	var div = copy_grid(w)
	for x in range(1, squares_qtd.y+1):
		for y in range(1, squares_qtd.x+1):
			#iterate over columns which is equivalent to iterate over x axis
			div[x][y].pressure = divergent_at_point(w, x, y)
	
	return div

func divergent_at_point(w, x, y):
	var dx = tile_size.x

	var right = w[x][y+1].velocity.x 
	var left = w[x][y-1].velocity.x 
	var up = w[x+1][y].velocity.y 
	var down = w[x-1][y].velocity.y
	
	return -0.5 * dx * (right - left) + (up - down)

func gradient(grid: Array):
	var grad = grid.duplicate(true)
	for x in range(1, squares_qtd.y+1):
		for y in range(1, squares_qtd.x+1):
			grad[x][y] = gradient_at_point(grid, x, y)
	
	return grad

func gradient_at_point(grid, x, y):
	var dx = tile_size.x
	var dy = tile_size.y 

	var right = grid[x][y+1]
	var left = grid[x][y-1]
	var up = grid[x+1][y]
	var down = grid[x-1][y]
	
	return 0.5 * Vector2((right-left)/dx, (up-down)/dy)

func project(w):
	#Refactor make func
	var x0 = []
	for x in range(squares_qtd.y+2):
		x0.append([])
		for _y in range(squares_qtd.x+2):
			x0[x].append(0)
	
	var q = poisson_solver(divergent(w), x0, 10e-5)
	var grad_q = gradient(q)
	for x in range(1, squares_qtd.y+1):
		for y in range(1, squares_qtd.x+1):
			w[x][y].velocity -= grad_q[x][y]
			w[x][y].check_vel(MAX_VELOCITY)
			w[x][y].pressure = rho/(0.5*(tile_size.x+tile_size.y)) * q[x][y]

func update_grid(new_field):
	for x in range(squares_qtd.y+2):
		for y in range(squares_qtd.x+2):
			grid_vectors[x][y].pressure = new_field[x][y].pressure
			grid_vectors[x][y].velocity = new_field[x][y].velocity

func update_field(delta):
	var w = copy_grid(grid_vectors)
	add_force(w, delta)
	update_boundary(w)
	advect(w, delta)
	update_boundary(w)
	diffuse(w)
	project(w)
	update_boundary(w)
	return w

func _process(delta):
	timer += delta
	
	if get_parent().interface_visible:
		return
	
	if Input.is_action_pressed("left_click"):
		if timer >= 0.1:
			add_particle()
			timer = 0
	
	var new_field = update_field(delta)
	update_grid(new_field)

func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	if show_vectors:
		$vector_visualizer.show()
	else:
		$vector_visualizer.hide()

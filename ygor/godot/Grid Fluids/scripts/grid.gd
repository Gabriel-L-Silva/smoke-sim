extends Node2D


var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(16, 9)
var tile_size 		= Vector2(grid_size.x/squares_qtd.x, grid_size.y/squares_qtd.y)
var show_vectors 	= false
var show_grid 		= false
var timer = 0
var vel_min = pow(2, 63)
var vel_max = -vel_min

var gravity = Vector2(0, 9.81)
var sub_steps = squares_qtd.x #random value, maybe be lowered for performance improvement

var Particle = preload("res://scenes/particle.tscn")
var Vector = preload("res://scenes/vector.tscn")

var grid_vectors = []
var particles = []

func get_velocity(pos):
	return Vector2(pos.x, -pos.y)

func get_pressure(pos):
	return pos.x+pos.y

func _ready():
	for x in range(squares_qtd.x):
		grid_vectors.append([])
		for y in range(squares_qtd.y):
			# middle point
			var pos = Vector2(x*tile_size.x+tile_size.x/2, y*tile_size.y+tile_size.y/2)
			var vector = Vector.instance()
			vector.position = pos
			vector.velocity = get_velocity(pos)
			
			if vector.velocity.length() > vel_max:
				vel_max = vector.velocity.length()
			
			if vector.velocity.length() < vel_min:
				vel_min = vector.velocity.length() 
			
			vector.pressure = get_pressure(pos) 
			grid_vectors[x].append(vector)
			$vector_visualizer.add_child(vector)
	
func bilinear_interpolation(pos) -> Vector2:
	var dx = tile_size.x
	var dy = tile_size.y
	
	var x0 = 0.5 * dx
	var y0 = 0.5 * dy
	
	var i = int((pos.x-x0)/dx)
	var j = int((pos.y-y0)/dy)
	
	if pos.x - x0 < 0:
		pos.x = x0
	
	if pos.y - y0 < 0:
		pos.y = y0
	
	if i+1 >= squares_qtd.x:
		pos.x = grid_size.x - x0
		i -= squares_qtd.x
	
	if j+1 >= squares_qtd.y:
		pos.y = grid_size.y - y0
		j -= squares_qtd.y

	var xi = grid_vectors[i][j].position.x
	var yj = grid_vectors[i][j].position.y
	
	var x1 = xi
	var y1 = yj
	var q11 = grid_vectors[i][j]
	
	var x2 = xi + dx
	var _y1 = yj
	var q21 = grid_vectors[i+1][j]
	
	var _x1 = xi
	var y2 = yj + dy
	var q12 = grid_vectors[i][j+1]
	
	var _x2 = xi + dx
	var _y2 = yj + dy
	var q22 = grid_vectors[i+1][j+1]

	return (q11.velocity * (x2 - pos.x) * (y2 - pos.y) +
		q21.velocity * (pos.x - x1) * (y2 - pos.y) +
		q12.velocity * (x2 - pos.x) * (pos.y - y1) +
		q22.velocity * (pos.x - x1) * (pos.y - y1)
		) / ((x2 - x1) * (y2 - y1))

func add_particle():
	var pos = get_global_mouse_position()
	var new_particle = Particle.instance()
	new_particle.position = pos
	particles.append(new_particle)
	add_child(new_particle)
	print('particle added')

func external_forces():
	return gravity #+ bouancy (+ mouse_reppelant_force)

func add_force(w0: Array, delta):
	var w1 = w0.duplicate(true)
	for x in range(squares_qtd.x):
		for y in range(squares_qtd.y):
			w1[x][y].velocity = w1[x][y].velocity + delta * external_forces()
	return w1

func advect(w1: Array, timestep):
	var w2 = w1.duplicate(true)
	var s = timestep/sub_steps
	for x in range(squares_qtd.x):
		for y in range(squares_qtd.y):
			var pos = grid_vectors[x][y].position
			for _s in range(sub_steps):
				pos = pos - s*bilinear_interpolation(pos)
			w2[x][y].velocity = bilinear_interpolation(pos)
	return w2

func diffuse(w2, delta):
	return w2 #Passthrough function because smoke has no viscosity

func poisson_solver(vec_field):
	#TODO implement poisson solver
	pass

func grad(vec_field):
	#TODO implement finite diferences grad (all cells or single cell? as is function project is using cell grad)
	pass
	
func project(w3, delta):
	var q = poisson_solver(w3)
	var w4 = w3.duplicate(true)
	for x in range(squares_qtd.x):
		for y in range(squares_qtd.y):
			w4[x][y].velocity = w3[x][y].velocity - grad(q)

func update_field(delta):
	var w0 = grid_vectors
	var w1 = add_force(w0, delta)
	var w2 = advect(w1, delta)
	var w3 = diffuse(w2, delta)
	#var w4 = project(w3, delta)
	return w3
func _process(delta):
	timer += delta
	
	grid_vectors = update_field(delta/5) #division so animation is slower
	#update vector drawing each frame GAMBITO
	#TODO improve vector drawing update method
	$vector_visualizer.hide()
	$vector_visualizer.show()
		
	if Input.is_action_pressed("left_click") and not get_parent().interface_visible:
		if timer >= 0.05:
			add_particle()
			timer = 0
	
	for p in particles:
		if p.must_remove:
			particles.erase(p)
			remove_child(p)
		else:
			p.vel = bilinear_interpolation(p.position)

func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	if show_vectors:
		$vector_visualizer.show()
	else:
		$vector_visualizer.hide()

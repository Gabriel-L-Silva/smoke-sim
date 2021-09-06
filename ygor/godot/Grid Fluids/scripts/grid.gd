extends Node2D

var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(32, 32)
var tile_size 		= Vector2(grid_size.x/squares_qtd.x, grid_size.y/squares_qtd.y)
var show_vectors 	= false
var show_grid 		= false
var timer = 0.0

var rho = 1.0
var gravity = Vector2(0, -9.81) 
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

func get_velocity(pos):
	return Vector2(pos.x/100, 0)

func get_pressure(pos):
	return 1.0

func copy_vector(obj):
	var vec = VectorClass.new()
	vec.pressure = obj.pressure
	vec.velocity = obj.velocity
	vec.pos = obj.pos
	return vec

func _ready():
	$native_lib.tile_size = tile_size
	$native_lib.grid_size = grid_size
	
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
	
	$native_lib.vector_size = Vector2(squares_qtd.y+2, squares_qtd.x+2)
	
	print($native_lib.grid_size)
	print($native_lib.tile_size)
	print($native_lib.vector_size)

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

func add_particle():
	var pos = get_global_mouse_position()
	var new_particle = Particle.instance()
	new_particle.position = pos
	particles.append(new_particle)
	add_child(new_particle)
	print('particle added')

func external_forces():
	return gravity # + bouancy (+ mouse_reppelant_force)

func update_grid_native(new_field):
	for x in range(squares_qtd.y+2):
		for y in range(squares_qtd.x+2):
			grid_vectors[x][y].pressure = new_field[x][y][0]
			grid_vectors[x][y].velocity = new_field[x][y][1]

func _process(delta):
	timer += delta
	
	if get_parent().interface_visible:
		return
	
	if Input.is_action_pressed("left_click"):
		if timer >= 0.1:
			add_particle()
			timer = 0
	
	var new_field = $native_lib.update_field(delta, grid_vectors, external_forces())
	update_grid_native(new_field)
	
func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	if show_vectors:
		$vector_visualizer.show()
	else:
		$vector_visualizer.hide()

extends Node2D

var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(64,64)
var tile_size 		= Vector2(grid_size.x/squares_qtd.x, grid_size.y/squares_qtd.y)
var show_vectors 	= false
var show_grid 		= false
var show_pressure 	= false
var show_density 	= false
var timer = 0.0
var mouse_inside = true
var minmax_vel = Vector2(0,1)

var rho = 1.0
var gravity = Vector2(0, 9.81) 
var sub_steps = 1 # random value, maybe be lowered for performance improvement
var MAX_VELOCITY = 500

var Particle = preload("res://scenes/particle.tscn")
var Vector = preload("res://scenes/vector.tscn")

var grid_vectors = []
var particles = []

class VectorClass:
	var pressure: float
	var density: float
	var velocity: Vector2
	var pos: Vector2

func get_velocity(_pos):
	return Vector2(-pow(_pos.y-grid_size.y/2.0,3)-9*(_pos.y-grid_size.y/2.0), -pow(_pos.x-grid_size.x/2.0,3)-9*(_pos.x-grid_size.x/2.0))/10000
#	return Vector2(cos(_pos.x+_pos.y), sin(_pos.x*_pos.y))*100
#	return Vector2(-_pos.y+grid_size.y/2, _pos.x-grid_size.x/2)
#   return _pos
#	return Vector2(10,10)

func get_pressure(_pos):
	return (_pos.x+_pos.y)/100 if _pos.x >= grid_size.x/2-2*tile_size.x and _pos.y >= grid_size.y/2-2*tile_size.y and _pos.x <= grid_size.x/2+2*tile_size.x and _pos.y <= grid_size.y/2+2*tile_size.y else 0
#	return (_pos.x/1280/2+_pos.y/720/2)

func get_density(_pos):
#	return 1 if _pos.x >= grid_size.x - tile_size.x else 0 
#	return 1 if _pos.x >= grid_size.x/2-2*tile_size.x and _pos.y >= grid_size.y/2-2*tile_size.y and _pos.x <= grid_size.x/2+2*tile_size.x and _pos.y <= grid_size.y/2+2*tile_size.y else 0.0
	return 1 if _pos.x >= grid_size.x/2-3*tile_size.x and _pos.y >= grid_size.y/2-3*tile_size.y and _pos.x <= grid_size.x/2+3*tile_size.x and _pos.y <= grid_size.y/2+3*tile_size.y else 0.0
#	return 1

func copy_vector(obj):
	var vec = VectorClass.new()
	vec.pressure = obj.pressure
	vec.density = obj.density
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
			var pos = Vector2((y+1) * tile_size.x + tile_size.x/2.0, (x+1) * tile_size.y + tile_size.y/2.0)
			var vector = Vector.instance()
			vector.pos = pos
			vector.velocity = get_velocity(pos)
			vector.pressure = get_pressure(pos)
			vector.density = get_density(pos)
			grid_vectors[x].append(vector)
			$vector_visualizer.add_child(vector)
		
		# vertical
		var front = copy_vector(grid_vectors[x][0])
		var back = copy_vector(grid_vectors[x][-1])
		front.pos.x = tile_size.x/2.0
		front.velocity *= -1 
		back.pos.x += tile_size.x
		back.velocity *= -1
		grid_vectors[x].push_front(front)
		grid_vectors[x].push_back(back)
	
	# horizontal
	var front_list = []
	var back_list = []
	for y in range(squares_qtd.x+2):
		var vec = copy_vector(grid_vectors[0][y])
		vec.pos.y = tile_size.y/2.0
		if y > 0:
			vec.velocity *= -1
		front_list.append(vec)
	for y in range(squares_qtd.x+2):
		var vec = copy_vector(grid_vectors[-1][y])
		vec.pos.y += tile_size.y
		if y>0:
			vec.velocity *= -1
		back_list.append(vec)
	grid_vectors.push_front(front_list)
	grid_vectors.push_back(back_list)
	
	$native_lib.vector_size = Vector2(squares_qtd.y+2, squares_qtd.x+2)
	minmax_vel = get_minmax_velocity()

func _notification(what):
	match what:
		MainLoop.NOTIFICATION_WM_MOUSE_EXIT:
#			print("Mouse has left game window")
			mouse_inside = false
		MainLoop.NOTIFICATION_WM_MOUSE_ENTER:
#			print("Mouse has entered the game window")
			mouse_inside = true

func get_minmax_velocity():
	var result = $native_lib.get_minmax_velocity(grid_vectors);
	return result
	
func get_minmax_pressure():
	var result = $native_lib.get_minmax_pressure(grid_vectors);
	return result

func bilinear_interpolation_vel(pos):
	var result = $native_lib.bilinear_interpolation_grid(grid_vectors, pos, false);
	return result

func bilinear_interpolation_press(pos):
	var result = $native_lib.bilinear_interpolation_grid(grid_vectors, pos, true);
	return result.x

func bilinear_interpolation_density(pos):
	var result = $native_lib.bilinear_interpolation_grid(grid_vectors, pos, true);
	return result.y

func add_particle(pos):
	var new_particle = Particle.instance()
	new_particle.position = pos
	particles.append(new_particle)
	add_child(new_particle)

func external_forces():
	return gravity # + bouancy (+ mouse_reppelant_force)

func _process(delta):
	timer += delta
#	OS.delay_msec(300)
	var pos = get_global_mouse_position()
	if Input.is_action_pressed("left_click"):
		if timer >= 0.05 and not get_parent().check_inter_col(pos):
			add_particle(pos)
			timer = 0
	
	if get_parent().interface_visible:
		return 
	
	$native_lib.mouse_pos = pos if mouse_inside else Vector2(-1,-1)
	var check = $native_lib.update_field(delta, grid_vectors, external_forces())
#	print(check)
	$native_lib.update_particles(grid_vectors, particles, delta)
	minmax_vel = get_minmax_velocity()

func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_pressure_signal():
	show_pressure = not show_pressure
	$pressure_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	if show_vectors:
		$vector_visualizer.show()
	else:
		$vector_visualizer.hide()

func _on_interface_show_density_signal():
	show_density = not show_density
	$density_visualizer.update()

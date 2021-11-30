extends Node2D

var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(32, 32)
var tile_size 		= Vector2(grid_size.x/squares_qtd.x, grid_size.y/squares_qtd.y)
var show_vectors 	= false
var show_grid 		= false
var timer = 0.0
var mouse_pos = Vector2(0, 0)
var use_bilinear_interp = false
var aditional_tiles: int

var rho = 1.0
var gravity = Vector2(0, -0.981) 
var sub_steps = 1 # one backtrack only usualy works best
var MAX_VELOCITY = 120

var Particle = preload("res://scenes/particle.tscn")
var Vector = preload("res://scenes/vector.tscn")

var grid_vectors = []
var particles = []

class VectorClass:
	var pressure: float
	var velocity: Vector2
	var pos: Vector2

func get_velocity(_pos):
#	return Vector2(_pos.x, _pos.y)/100
	return Vector2(0,0)

func get_pressure(_pos):
	return 10 #(_pos.x+_pos.y)/100

func copy_vector(obj):
	var vec = VectorClass.new()
	vec.pressure = obj.pressure
	vec.velocity = obj.velocity
	vec.pos = obj.pos
	return vec

func _ready():
	$native_lib.tile_size = tile_size
	$native_lib.grid_size = grid_size
	$native_lib.max_speed = MAX_VELOCITY
	$native_lib.rho_const = rho
	$native_lib.bilinear_interp = use_bilinear_interp
	
	if use_bilinear_interp: aditional_tiles=2
	else: aditional_tiles=4
	
	for x in range(squares_qtd.y):
		grid_vectors.append([])
		for y in range(squares_qtd.x):
			# middle point
			var pos: Vector2
			if use_bilinear_interp:
				pos = Vector2((y+1) * tile_size.x + tile_size.x/2.0, (x+1) * tile_size.y + tile_size.y/2.0)
			else:
				pos = Vector2((y+2) * tile_size.x + tile_size.x/2.0, (x+2) * tile_size.y + tile_size.y/2.0)
			
			var vector = Vector.instance()
			vector.pos = pos
			vector.velocity = get_velocity(pos)
			vector.pressure = get_pressure(pos) 
			grid_vectors[x].append(vector)
			$vector_visualizer.add_child(vector)
			
		# vertical
		if use_bilinear_interp:
			var front = copy_vector(grid_vectors[x][0])
			var back = copy_vector(grid_vectors[x][-1])
			front.pos.x = tile_size.x/2.0
			front.velocity *= -1 
			back.pos.x += tile_size.x
			back.velocity *= -1
			grid_vectors[x].push_front(front)
			grid_vectors[x].push_back(back)
		else:
			# if camull need second layer
			var front1 = copy_vector(grid_vectors[x][0])
			var back1 = copy_vector(grid_vectors[x][-1])
			var front2 = copy_vector(grid_vectors[x][0])
			var back2 = copy_vector(grid_vectors[x][-1])
			front1.pos.x = tile_size.x/2.0 + tile_size.x
			front1.velocity *= -1 
			back1.pos.x += tile_size.x
			back1.velocity *= -1
			front2.pos.x = tile_size.x/2.0
			front2.velocity *= -1 
			back2.pos.x += tile_size.x + tile_size.x
			back2.velocity *= -1
			grid_vectors[x].push_front(front1)
			grid_vectors[x].push_back(back1)
			grid_vectors[x].push_front(front2)
			grid_vectors[x].push_back(back2)
	
	# horizontal
	if use_bilinear_interp:
		var front_list = []
		var back_list = []
		for y in range(squares_qtd.x + aditional_tiles):
			var vec = copy_vector(grid_vectors[0][y])
			vec.pos.y = tile_size.y/2.0
			if y > 0:
				vec.velocity *= -1
			front_list.append(vec)
		for y in range(squares_qtd.x + aditional_tiles):
			var vec = copy_vector(grid_vectors[-1][y])
			vec.pos.y += tile_size.y
			if y>0:
				vec.velocity *= -1
			back_list.append(vec)
		grid_vectors.push_front(front_list)
		grid_vectors.push_back(back_list)
	else:
		# if camull need second layer
		var front_list1 = []
		var back_list1 = []
		var front_list2 = []
		var back_list2 = []
		for y in range(squares_qtd.x + aditional_tiles):
			var vec1 = copy_vector(grid_vectors[0][y])
			var vec2 = copy_vector(grid_vectors[0][y])
			vec1.pos.y = tile_size.y/2.0 + tile_size.y
			vec2.pos.y = tile_size.y/2.0
			if y > 0:
				vec1.velocity *= -1
				vec2.velocity *= -1
			front_list1.append(vec1)
			front_list2.append(vec2)
		for y in range(squares_qtd.x + aditional_tiles):
			var vec1 = copy_vector(grid_vectors[-1][y])
			var vec2 = copy_vector(grid_vectors[-1][y])
			vec1.pos.y += tile_size.y
			vec2.pos.y += tile_size.y + tile_size.y
			if y>0:
				vec1.velocity *= -1
				vec2.velocity *= -1
			back_list1.append(vec1)
			back_list2.append(vec2)
		grid_vectors.push_front(front_list1)
		grid_vectors.push_back(back_list1)
		grid_vectors.push_front(front_list2)
		grid_vectors.push_back(back_list2)
	
	$native_lib.vector_size = Vector2(squares_qtd.y+aditional_tiles, squares_qtd.x+aditional_tiles)

func add_particle(pos):
	if pos.x < 0 or pos.x > grid_size.x or pos.y < 0 or pos.y > grid_size.y:
		return
		
	var new_particle = Particle.instance()
	new_particle.position = pos
	particles.append(new_particle)
	add_child(new_particle)

func external_forces():
#	var force = gravity
	var force = Vector2(0, 0)
	
	if Input.is_action_pressed("ui_left"):
		force += Vector2(-50, 0)
	if Input.is_action_pressed("ui_right"):
		force += Vector2(50, 0)
	if Input.is_action_pressed("ui_up"):
		force += Vector2(0, -50)
	if Input.is_action_pressed("ui_down"):
		force += Vector2(0, 50)
	
	return force
	
func add_mouse_repel(pos, prev, delta):
	if pos.x < 0 or pos.x > grid_size.x or pos.y < 0 or pos.y > grid_size.y or (
		prev.x < 0 or prev.x > grid_size.x or prev.y < 0 or prev.y > grid_size.y):
		return
	
	var i = floor( pos.y / tile_size.y) + 1;
	var j = floor( pos.x / tile_size.x) + 1;
	
	if i == squares_qtd.x:
		i -= 1
	if j == squares_qtd.y:
		j -= 1
	
	var v: Vector2 = pos-prev
	v = v.normalized()
	
	var add = MAX_VELOCITY
	if not use_bilinear_interp:
		add /= 4
		
	grid_vectors[i][j].velocity += v * add
	grid_vectors[i+1][j].velocity += v * add
	grid_vectors[i][j+1].velocity += v * add
	grid_vectors[i+1][j+1].velocity += v * add


func _process(delta):
	timer += delta
	var cur_mouse = get_global_mouse_position()
	
	if Input.is_action_pressed("left_click"):
		if timer >= 0.05 and not get_parent().check_inter_col(cur_mouse):
#			print(floor((cur_mouse.x - tile_size.x/2.0) / tile_size.x) + 2)
			add_particle(cur_mouse)
			timer = 0

	if Input.is_action_pressed("right_click"):
		add_mouse_repel(cur_mouse, mouse_pos, delta)
	
	mouse_pos = cur_mouse
	
	if get_parent().interface_visible:
		return
	
	$native_lib.update_field(delta, grid_vectors, particles, external_forces())


func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	if show_vectors:
		$vector_visualizer.show()
	else:
		$vector_visualizer.hide()

extends Node2D


var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(16, 9)
var tile_size 		= Vector2(grid_size.x/squares_qtd.x, grid_size.y/squares_qtd.y)
var show_vectors 	= false
var show_grid 		= false
var timer = 0

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
			vector.pressure = get_pressure(pos) 
			grid_vectors[x].append(vector)
			$vector_visualizer.add_child(vector)

func bilinear_interpolation(pos):
	#if(pos.x < 0.5 * tile_size.x or pos.y < 0.5 * tile_size.y
	#	or pos.x > grid_size.x - 0.5 * tile_size.x
	#	or pos.y > grid_size.y - 0.5 * tile_size.y):
	#	print("no border case")
	#	return

	var dx = tile_size.x
	var dy = tile_size.y
	
	var x0 = 0.5 * dx
	var y0 = 0.5 * dy
	
	var i = int((pos.x-x0)/dx)
	var j = int((pos.y-y0)/dy)

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

	if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
		print('points do not form a rectangle')
		return
	#if not (x1 <= pos.x and pos.x <= x2) or not (y1 <= pos.y and pos.y <= y2):
	#	print('(x, y) not within the rectangle')
	#	return

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
	
	print(pos)
	print(bilinear_interpolation(pos))

func _process(delta):
	timer += delta
	if Input.is_action_pressed("left_click") and not get_parent().interface_visible:
		if timer >= 0.15:
			add_particle()
			timer = 0
	
	for p in particles:
		p.vel = get_velocity(p.position)

func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	if show_vectors:
		$vector_visualizer.show()
	else:
		$vector_visualizer.hide()

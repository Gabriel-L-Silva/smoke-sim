extends TileMap


var tile_size 		= Vector2(32, 32)
var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(round(grid_size.x/tile_size.x), round(grid_size.y/tile_size.y))
var show_vectors 	= false
var show_grid 		= false
#var img_grid = Image.new()

var grid_vectors = []
var particles = []

class Vect:
	var velocity: Vector2
	var pos: Vector2
	var pressure: float
	
	func _init(velocity_, pos_, pressure_):
		velocity = velocity_
		pos = pos_
		pressure = pressure_
	
		
class Particle:
	var p: float
	var vel: Vector2
	var pos: Vector2
		
	func _init(pos_, velocity_, pressure_):
		vel = velocity_
		pos = pos_
		p = pressure_


func get_velocity(pos):
	return Vector2(0, pos.y)

func get_pressure(pos):
	return pos.x+pos.y

func _ready():
	#img_grid.create(squares_qtd.x, squares_qtd.y, false, Image.FORMAT_RGBA8)
	#img_grid.lock()
	
	for x in range(squares_qtd.x):
		grid_vectors.append([])
		for y in range(squares_qtd.y):
			# middle point
			var pos = Vector2(x*tile_size.x+tile_size.x/2, y*tile_size.y+tile_size.y/2)
			var velocity = get_velocity(pos)
			var p = get_pressure(pos) 
			grid_vectors[x].append(Vect.new(velocity, pos, p))
			#img_grid.set_pixel(x, y, Color(vel_x, vel_y, p, 1))
			
	#img_grid.unlock()
	#img_grid.resize(1280, 720, Image.INTERPOLATE_BILINEAR)
	#img_grid.save_png("res://textures/grid_info.png")
	#$smoke.set_material(load("res://shaders/grid.tres"))
	#var t = ImageTexture.new()
	#t.create_from_image(img_grid, 0)
	#$smoke.get_material().set_shader_param("noise_img", t)

func interpolate(pos):
	if(pos.x < 0.5 * tile_size.x or pos.y < 0.5 * tile_size.y
		or pos.x > grid_size.x - 0.5 * tile_size.x
		or pos.y > grid_size.y - 0.5 * tile_size.y):
		print("no border case")
		return

	var dx = tile_size.x
	var dy = tile_size.y
	
	var x0 = 0.5 * dx
	var y0 = 0.5 * dy
	var i = int((pos.x-x0)/dx)
	var j = int((pos.y-y0)/dy)

	var xi = grid_vectors[i][j].pos.x
	var yj = grid_vectors[i][j].pos.y
	
	var x1 = xi
	var y1 = yj
	var q11 = grid_vectors[j][i]
	
	var x2 = xi + dx
	var _y1 = yj
	var q21 = grid_vectors[j][i+1]
	
	var _x1 = xi
	var y2 = yj + dy
	var q12 = grid_vectors[j+1][i]
	
	var _x2 = xi + dx
	var _y2 = yj + dy
	var q22 = grid_vectors[j+1][i+1]

	if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
		print('points do not form a rectangle')
		return
	if not (x1 <= pos.x and pos.x <= x2) or not (y1 <= pos.y and pos.y <= y2):
		print('(x, y) not within the rectangle')
		return

	return (q11.velocity * (x2 - pos.x) * (y2 - pos.y) +
		q21.velocity * (pos.x - x1) * (y2 - pos.y) +
		q12.velocity * (x2 - pos.x) * (pos.y - y1) +
		q22.velocity * (pos.x - x1) * (pos.y - y1)
		) / ((x2 - x1) * (y2 - y1))

func _process(delta):
	if Input.is_action_just_pressed("left_click") and not get_parent().interface_visible:
		var pos = get_global_mouse_position()
		particles.append(Particle.new(pos, Vector2(0,0), 0))
		var v = interpolate(pos)
		print(v)
		print(pos)
		
	for p in particles:
		pass

func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	$vector_visualizer.update()


extends TileMap


var tile_size 		= get_cell_size()
var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(round(grid_size.x/tile_size.x), round(grid_size.y/tile_size.y))
var show_vectors 	= false
var show_grid 		= false
var img_grid = Image.new()

var grid_vectors = []


class Vect:
	var velocity: Vector2
	var pos: Vector2
	var pressure: float
	
	func _init(velocity_, pressure_, pos_):
		velocity 	= velocity_
		pressure 	= pressure_
		pos 	 	= pos_


func _ready():
	img_grid.create(squares_qtd.x, squares_qtd.y, false, Image.FORMAT_RGBA8)
	img_grid.lock()
	
	for x in range(squares_qtd.x):
		grid_vectors.append([])
		for y in range(squares_qtd.y):
			# middle point
			var pos = Vector2(x*tile_size.x+tile_size.x/2, y*tile_size.y+tile_size.y/2)
			var vel_x = (randi()%255+1) # * -1
			var vel_y = (randi()%255+1) # * -1
			var p 	  = (randi()%255+1)
			grid_vectors[x].append(Vect.new(Vector2(vel_x, vel_y), p, pos))
			img_grid.set_pixel(x, y, Color(vel_x, vel_y, p, 1))
			
	img_grid.unlock()
	#img_grid.resize(1280, 720, Image.INTERPOLATE_BILINEAR)
	img_grid.save_png("res://textures/grid_info.png")
	$smoke.set_material(load("res://shaders/grid.tres"))
	var t = ImageTexture.new()
	t.create_from_image(img_grid, 0)
	$smoke.get_material().set_shader_param("noise_img", t)

func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	$vector_visualizer.update()


extends Node2D


onready var grid = get_parent()
var dynamic_font = DynamicFont.new()

func _ready():
	dynamic_font.font_data = load("res://textures/04B_19__.ttf")
	dynamic_font.size = 16

func _process(delta):
	update()

func normalize(value, minmax):
	return value
#	return (value - minmax.x)/(minmax.y - minmax.x)
	
func _draw():
	if not grid.show_grid:
		return
		
#	# You'll have to get thoose the way you want
#	var array_width = grid.squares_qtd.x+1
#	var array_height = grid.squares_qtd.y+1
#
#	var byte_array = PoolByteArray()
#	for x in range(1,grid.squares_qtd.y+2):
#		var line = []
#		for y in range(1,grid.squares_qtd.x+2):
#			var pos = grid.grid_vectors[x][y].pos - 0.5*grid.tile_size
#			var position = pos-grid.tile_size
#			line.append(grid.bilinear_interpolation_press(position))
##			draw_circle(position, 1, Color(1,0,0,1))
#		byte_array.append_array(PoolByteArray(line))
#
#	# The following is used to convert the array into a Texture
#	var img = Image.new()
#
#	# I don't want any mipmaps generated : use_mipmaps = false
#	img.create_from_data(array_width, array_height, false, Image.FORMAT_R8, byte_array)
#
#	var texture = ImageTexture.new()
#
#	# Override the default flag with 0 since I don't want texture repeat/filtering/mipmaps/etc
#	texture.create_from_image(img, 0)
#
#	# Upload the texture to my shader
#	grid.material.set_shader_param("my_array", texture)
	var LINE_COLOR = Color(0, 0, 0)
	var LINE_WIDTH = 2

	var limit_x = grid.squares_qtd.x * grid.tile_size.x
	var limit_y = grid.squares_qtd.y * grid.tile_size.y

	var vertex_pressure = []
	for x in range(1,grid.squares_qtd.y+2):
		vertex_pressure.append([])
		for y in range(1,grid.squares_qtd.x+2):

			var pos = grid.grid_vectors[x][y].pos - 0.5*grid.tile_size
			var position = pos-grid.tile_size
			var vec = grid.VectorClass.new()
			vec.pressure = grid.bilinear_interpolation_press(position)
#			if(x==5 and y==5): 
#				vec.pressure = 255
			vec.pos = position
			vertex_pressure[x-1].append(vec)
#			draw_circle(position + 0.5*grid.tile_size, 1, Color(1,1,1,vec.pressure/255.0))
	
	var minmax_presure = grid.get_minmax_pressure()
	#da pra otimizar, esta calculando o mesmo ponto várias vezes
	for x in range(vertex_pressure.size()-1):
		for y in range(vertex_pressure[0].size()-1):
			var position = vertex_pressure[x][y].pos
			var lb = position
			var lt = lb + Vector2(0,grid.tile_size.y)
			var rb = lb + Vector2(grid.tile_size.x, 0)
			var rt = rb + Vector2(0,grid.tile_size.y)
			var verts = PoolVector2Array([lb, rb, rt, lt])
			var colors = PoolColorArray([Color(1,1,1, normalize(vertex_pressure[x][y].pressure, minmax_presure)),
										 Color(1,1,1, normalize(vertex_pressure[x+1][y].pressure, minmax_presure)),
										 Color(1,1,1, normalize(vertex_pressure[x+1][y+1].pressure, minmax_presure)),
										 Color(1,1,1, normalize(vertex_pressure[x][y+1].pressure, minmax_presure))])
##			var colors = PoolColorArray()
#			var uvs = PoolVector2Array()
#			for vert in verts:
##				colors.append(Color(1,1,1, vertex_pressure[x][y].pressure/255.0))
#				uvs.append(Vector2(position.x/grid.grid_size.x, position.y/grid.grid_size.y))
			draw_primitive(verts, colors, verts)
#			draw_string(dynamic_font, lb+Vector2(0,dynamic_font.size), "%.2f" %vertex_pressure[x][y].pressure,Color(1,1,0))
#			if(x==4 and y==4):
#				draw_string(dynamic_font, rb, "%.2f" %vertex_pressure[x+1][y].pressure,Color(1,1,0))
#				draw_string(dynamic_font, rt, "%.2f" %vertex_pressure[x+1][y+1].pressure,Color(1,1,0))
#				draw_string(dynamic_font, lt, "%.2f" %vertex_pressure[x][y+1].pressure,Color(1,1,0))
#				draw_circle(lb,3,Color(1,0,0))
#				draw_circle(rb,3,Color(0,1,0))
#				draw_circle(rt,3,Color(0,0,1))
#				draw_circle(lt,3,Color(1,0,1))
#				print(' ')
#			draw_rect(Rect2(position,grid.tile_size), Color(1,1,1,1), false, 1.0, true)
#			var colors = PoolColorArray([Color(1,1,1,1),Color(1,0,1,1),Color(1,1,0,1),Color(1,0,0,1)])
#			var colors = PoolColorArray()
#			var uvs = PoolVector2Array()
#			for vert in verts:
#				colors.append(Color(1,1,1, grid.bilinear_interpolation_press(pos)/255.0))
#				uvs.append(Vector2(1,1))
#			draw_primitive(verts, colors, uvs)
#			draw_rect(Rect2(position-grid.tile_size/2.0,grid.tile_size), Color(1,1,1,1))
#			draw_circle(position,1,Color(1,0,0,1))
			
#	# cols
#	for x in range(1,grid.squares_qtd.x+2):
#		var col_pos = x * grid.tile_size.x
#
#		draw_line(Vector2(col_pos, 0), 
#				  Vector2(col_pos, limit_y), LINE_COLOR, LINE_WIDTH, true)
#
#	# rows
#	for y in range(1,grid.squares_qtd.y+2):
#		var row_pos = y * grid.tile_size.y + LINE_WIDTH
#		draw_line(Vector2(0, row_pos), 
#				  Vector2(limit_x, row_pos), LINE_COLOR, LINE_WIDTH, true)


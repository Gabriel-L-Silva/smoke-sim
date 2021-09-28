extends Node2D

onready var grid = get_parent()
var dynamic_font = DynamicFont.new()

func _ready():
	dynamic_font.font_data = load("res://textures/04B_19__.ttf")
	dynamic_font.size = 12

func _process(delta):
	update()

func _draw():
	if not grid.show_density:
		return
#	var vertex = grid.get_density_primitive_vertex()
#	var colors = grid.get_density_primitive_colors()
	for x in range(0,grid.squares_qtd.y+1):
		for y in range(0,grid.squares_qtd.x+1):
			var position = grid.grid_vectors[x][y].pos
			var lb = position - grid.tile_size
			var lt = lb + Vector2(0,grid.tile_size.y)
			var rb = lb + Vector2(grid.tile_size.x, 0)
			var rt = rb + Vector2(0,grid.tile_size.y)
			var verts = PoolVector2Array([lb, rb, rt, lt])
			var colors = PoolColorArray([Color(1,1,1, grid.grid_vectors[x][y].density),
										 Color(1,1,1, grid.grid_vectors[x][y+1].density),
										 Color(1,1,1, grid.grid_vectors[x+1][y+1].density),
										 Color(1,1,1, grid.grid_vectors[x+1][y].density)])

			draw_primitive(verts, colors, verts)
#			if x==5 and y==5:
#				print('')
#			draw_string(dynamic_font, lb+Vector2(0,dynamic_font.size), "%.4f" %grid.grid_vectors[x][y].density,Color(1,0,0))
#			draw_string(dynamic_font, lb+Vector2(0,dynamic_font.size), "(%d,%d)" %[x,y] ,Color(0,1,0))
#			draw_circle(position, 1, Color(1,0,0,1))

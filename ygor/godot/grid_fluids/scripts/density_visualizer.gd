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
	var prim = grid.get_density_primitive()
	for x in range(0,prim.size(),2):
		var vertex = prim[x]
		var colors = prim[x+1]
		draw_primitive(vertex, colors, vertex)
#			if x==5 and y==5:
#				print('')
#			draw_string(dynamic_font, lb+Vector2(0,dynamic_font.size), "%.4f" %grid.grid_vectors[x][y].density,Color(1,0,0))
#			draw_string(dynamic_font, lb+Vector2(0,dynamic_font.size), "(%d,%d)" %[x,y] ,Color(0,1,0))
#			draw_circle(position, 1, Color(1,0,0,1))

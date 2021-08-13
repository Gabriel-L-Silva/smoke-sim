extends Node2D


onready var grid = get_parent()
var size = 10
var color = Color(255,255,255)

func _draw():
	if not grid.show_vectors:
		return 

	for x in grid.grid_vectors:
		for v in x:
			draw_vector(v)

func draw_vector(vect):
	var thick = ((abs(vect.velocity.x) + abs(vect.velocity.y)) / 2) / 255 * 2
	draw_line(vect.pos, get_end_point(vect), color, thick, true)
	draw_triangle(vect, thick)

func draw_triangle(vect, thick):
		var end_point = get_end_point(vect)
		var dir = vect.pos.direction_to(end_point)
		
		var a = end_point + dir * thick*2
		var b = end_point + dir.rotated(2*PI/3) * thick*2
		var c = end_point + dir.rotated(4*PI/3) * thick*2
		
		var points = PoolVector2Array([a, b, c])
		draw_polygon(points, PoolColorArray([color]), PoolVector2Array(), null, null, true)

func get_end_point(vect):
	return vect.pos + Vector2(size * (-1 if vect.velocity.x < 0 else 1) * (1 if vect.velocity.x else 0), 
							  size * (-1 if vect.velocity.y < 0 else 1) * (1 if vect.velocity.y else 0))


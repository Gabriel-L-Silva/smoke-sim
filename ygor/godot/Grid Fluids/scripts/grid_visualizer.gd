extends Node2D


onready var grid = get_parent()

func _process(delta):
	update()

func _draw():
	if not grid.show_grid:
		return

	var LINE_COLOR = Color(0, 0, 0)
	var LINE_WIDTH = 2

	var limit_x = grid.squares_qtd.x * grid.tile_size.x
	var limit_y = grid.squares_qtd.y * grid.tile_size.y

	#da pra otimizar, esta calculando o mesmo ponto várias vezes
	for x in range(1,grid.squares_qtd.y+1):
		for y in range(1,grid.squares_qtd.x+1):
			var pos = grid.grid_vectors[x][y].pos
			var position = pos-grid.tile_size
			if pos.x < 0: pos.x = 0
			if pos.y < 0: pos.y = 0 
			if pos.x > grid.grid_size.x: pos.x = grid.grid_size.x 
			if pos.y > grid.grid_size.y: pos.y = grid.grid_size.y
			var lb = position-grid.tile_size/2
			var lt = lb + Vector2(0,grid.tile_size.y)
			var rb = lb + Vector2(grid.tile_size.x, 0)
			var rt = rb + Vector2(0,grid.tile_size.y)
			var verts = PoolVector2Array([lb, rb, rt, lt])
#			var colors = PoolColorArray([Color(1,1,1,1),Color(1,0,1,1),Color(1,1,0,1),Color(1,0,0,1)])
			var colors = PoolColorArray()
			var uvs = PoolVector2Array()
			for vert in verts:
				colors.append(Color(1,1,1,grid.bilinear_interpolation_press(pos)/255.0))
				uvs.append(Vector2(1,1))
			draw_primitive(verts, colors, uvs)
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


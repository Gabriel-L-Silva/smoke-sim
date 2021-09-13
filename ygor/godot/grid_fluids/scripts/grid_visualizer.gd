extends Node2D

onready var grid = get_parent()

func _process(delta):
	update()

func _draw():
	if grid.show_grid:
		var LINE_COLOR = Color(1, 1, 1)
		var LINE_WIDTH = 2

		var limit_x = grid.squares_qtd.x * grid.tile_size.x
		var limit_y = grid.squares_qtd.y * grid.tile_size.y
		
		# cols
		for x in range(1,grid.squares_qtd.x+2):
			var col_pos = x * grid.tile_size.x

			draw_line(Vector2(col_pos, 0), 
					  Vector2(col_pos, limit_y), LINE_COLOR, LINE_WIDTH, true)

		# rows
		for y in range(1,grid.squares_qtd.y+2):
			var row_pos = y * grid.tile_size.y + LINE_WIDTH
			draw_line(Vector2(0, row_pos), 
					  Vector2(limit_x, row_pos), LINE_COLOR, LINE_WIDTH, true)

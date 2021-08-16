extends Node2D


onready var grid = get_parent()

func _draw():
	if not grid.show_grid:
		return

	var LINE_COLOR = Color(255, 255, 255)
	var LINE_WIDTH = 2

	var limit_x = grid.squares_qtd.x * grid.tile_size.x
	var limit_y = grid.squares_qtd.y * grid.tile_size.y

	# cols
	for x in range(grid.squares_qtd.x):
		var col_pos = x * grid.tile_size.x
		draw_line(Vector2(col_pos, 0), 
				  Vector2(col_pos, limit_y), LINE_COLOR, LINE_WIDTH, true)

	# rows
	for y in range(grid.squares_qtd.y):
		var row_pos = y * grid.tile_size.y + LINE_WIDTH
		draw_line(Vector2(0, row_pos), 
				  Vector2(limit_x, row_pos), LINE_COLOR, LINE_WIDTH, true)


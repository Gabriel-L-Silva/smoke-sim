extends Node2D


var velocity: Vector2
var pressure: float
var density: float
var pos: Vector2

onready var grid = get_parent().get_parent()
var dynamic_font = DynamicFont.new()

func _ready():
	dynamic_font.font_data = load("res://textures/04B_19__.ttf")
	dynamic_font.size = 12
	position.x = pos.x - grid.tile_size.x
	position.y = pos.y - grid.tile_size.y
	set_physics_process(false)

# warning-ignore:unused_argument
func _process(delta):
	if grid.show_vectors:
		update()

func normalize(value):
	# from 0.5 to 1 source: https://stats.stackexchange.com/a/281164
	return (value - grid.minmax_vel.x)/(grid.minmax_vel.y - grid.minmax_vel.x) * 0.5 + 0.5 if grid.minmax_vel.x != grid.minmax_vel.y else 1

func _draw():
	var l = normalize(velocity.length())
	var c = Color.from_hsv(170/360.0, l, 1)
	
#	if pos.x - grid.tile_size.x < 0 and pos.y - grid.tile_size.y < 0:
#		print('aq')
#		draw_string(dynamic_font,Vector2(grid.tile_size.x/2,grid.tile_size.y/2 + dynamic_font.size	),String(velocity),Color(1,0,0,1))
#	elif pos.x - grid.tile_size.x < 0:
#		print('aq')
#		draw_string(dynamic_font,Vector2(-grid.tile_size.x/2,pos.y-1.5*grid.tile_size.y),String(velocity),Color(1,0,0,1))
#	elif pos.y - grid.tile_size.y < 0:
#		print('aq')
#		draw_string(dynamic_font,Vector2(pos.x-1.5*grid.tile_size.x,grid.tile_size.y/2 + dynamic_font.size),String(velocity),Color(1,0,0,1))
#	else:
#	draw_string(dynamic_font,Vector2(0,0),String(velocity),Color(1,0,0,1))
	$Sprite.modulate = c
	$Sprite.rotation = velocity.angle() + PI/2
	$Sprite.scale = grid.tile_size/$Sprite.texture.get_size().length()
	$Sprite.scale *= l

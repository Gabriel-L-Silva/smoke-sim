extends Node2D


var velocity: Vector2
var pressure: float
var pos: Vector2

onready var grid = get_parent().get_parent()

func _ready():
	position.x = pos.x - grid.tile_size.x
	position.y = pos.y - grid.tile_size.y
	set_physics_process(false)

# warning-ignore:unused_argument
func _process(delta):
	update()

func normalize(value):
	# from 0.5 to 1 source: https://stats.stackexchange.com/a/281164
	return (value - grid.minmax_vel.x)/(grid.minmax_vel.y - grid.minmax_vel.x) * 0.5 + 0.5

func _draw():
	var l = normalize(velocity.length())
	var c = Color.from_hsv(170/360.0, l, 1)

	$Sprite.modulate = c
	$Sprite.rotation = velocity.angle() + PI/2
	$Sprite.scale = grid.tile_size/$Sprite.texture.get_size().length()
	$Sprite.scale *= l

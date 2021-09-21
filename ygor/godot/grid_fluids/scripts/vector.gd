extends Node2D


var velocity: Vector2
var pressure: float
var pos: Vector2

onready var grid = get_parent().get_parent()

func _ready():
	position = pos - grid.tile_size
	set_physics_process(false)

func _draw():
	var l = velocity.length()
	var c = Color.from_hsv(170/360.0, l/grid.MAX_VELOCITY, 1)

	$Sprite.modulate = c
	$Sprite.rotation = velocity.angle() + PI/2
	$Sprite.scale = grid.tile_size/$Sprite.texture.get_size().length()
